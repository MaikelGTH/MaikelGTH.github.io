/* Modelo 0D de ciclo cerrado: MEP / MEC / HCCI
 * Port en JavaScript de motor0d (Python/MATLAB). Mismas ecuaciones; el
 * integrador es RK4 de paso fijo (0,25°) en lugar de un paso adaptativo.
 *
 *   dp/dθ = (γ-1)/V · (dQc/dθ - dQw/dθ) - γ · p/V · dV/dθ,   γ = γ(T)
 */
(function (root) {
  'use strict';

  const AMB = { p: 1.013e5, R: 287, pEsc: 1.05 * 1.013e5 };

  const GEO = (() => {
    const B = 0.086, S = 0.086, L = 0.145, nCil = 4;
    const Ap = Math.PI * B * B / 4;
    return { B, S, L, nCil, Ap, Vd: Ap * S };
  })();

  // Parámetros base de cada motor (misma geometría para los tres)
  const MOTORES = {
    MEP: {
      nombre: 'MEP', rc: 10.5, g: 1.32, Tivc: 320, LHV: 43.0e6, AFR: 14.7, etaComb: 0.97,
      ctrl: 'cantidad', llama: 'premezcla', lambdaPC: 1.0,
      comb: { tipo: 'wiebe', soc: -22, dur: 50, a: 5, m: 2 },
      nMax: 6500, nEtaV: 4500, etaVmax: 0.92, xO2: 0.21
    },
    MEC: {
      nombre: 'MEC', rc: 17.5, g: 1.35, Tivc: 330, LHV: 42.6e6, AFR: 14.5, etaComb: 0.98,
      ctrl: 'calidad', llama: 'difusion', lambdaPC: 1.35,              // límite de humos
      comb: { tipo: 'doble', soc: -6, dur: 8, a: 5, m: 2, beta: 0.15, dur2: 55, m2: 0.9 },
      nMax: 4500, nEtaV: 2200, etaVmax: 0.90, xO2: 0.21
    },
    HCCI: {
      nombre: 'HCCI', rc: 15, g: 1.36, Tivc: 440, LHV: 43.0e6, AFR: 14.7, etaComb: 0.95,
      ctrl: 'hcci', llama: 'premezcla', lambdaPC: 1.0, lambdaMin: 2.0,
      comb: { tipo: 'wiebe', soc: null, dur: 12, a: 5, m: 2 },          // soc por Livengood-Wu
      combSI: { tipo: 'wiebe', soc: -8, dur: 55, a: 5, m: 2 },          // modo chispa
      TivcSI: 330, nMax: 6500, nEtaV: 4000, etaVmax: 0.88,
      xO2: 0.19,                                                        // dilución por residuales
      ca50Obj: 6                                                        // null = lazo abierto
    }
  };

  const OPC = { woschni: true, Tpared: 450, dth: 0.5, sub: 2, finLW: 30 };

  // ------------------------------------------------------------------ geometría
  function volumen(geo, rc, th) {
    const t = th * Math.PI / 180, a = geo.S / 2, R = geo.L / a;
    const Vc = geo.Vd / (rc - 1);
    const s = Math.sin(t), c = Math.cos(t);
    const raiz = Math.sqrt(R * R - s * s);
    return [Vc + geo.Ap * a * (R + 1 - c - raiz),
            geo.Ap * a * s * (1 + c / raiz) * Math.PI / 180];       // [m3], [m3/deg]
  }

  // ------------------------------------------------------- liberación de calor
  function wiebe(th, soc, dur, a, m) {
    const z = Math.max((th - soc) / dur, 0);
    const e = Math.exp(-a * Math.pow(z, m + 1));
    return [1 - e, a * (m + 1) / dur * Math.pow(z, m) * e];        // xb, dxb [1/deg]
  }

  function quemado(th, cb) {
    if (!cb) return [0, 0];
    if (cb.tipo === 'wiebe') return wiebe(th, cb.soc, cb.dur, cb.a, cb.m);
    const [x1, d1] = wiebe(th, cb.soc, cb.dur, cb.a, cb.m);         // Watson: premezcla
    const [x2, d2] = wiebe(th, cb.soc, cb.dur2, cb.a, cb.m2);       //         + difusión
    return [cb.beta * x1 + (1 - cb.beta) * x2, cb.beta * d1 + (1 - cb.beta) * d2];
  }

  const gammaT = T => 1.392 - 8.13e-5 * T;

  // ------------------------------------------------- estrategia de carga
  function estrategia(mot, carga) {
    let cb = mot.comb;
    if (mot.ctrl === 'cantidad') return { pIn: carga * AMB.p, lam: 1, cb, modo: 'Chispa' };
    if (mot.ctrl === 'calidad') {
      cb = Object.assign({}, cb, { beta: 0.10 + 0.35 * (1 - carga), dur2: 30 + 30 * carga });
      return { pIn: AMB.p, lam: mot.lambdaPC / carga, cb, modo: 'Difusión' };
    }
    if (carga <= 1 / mot.lambdaMin + 1e-9) return { pIn: AMB.p, lam: 1 / carga, cb, modo: 'HCCI' };
    return { pIn: carga * AMB.p, lam: 1, cb: mot.combSI, modo: 'Chispa' };
  }

  // Curva SUPUESTA (no calculada): parábola centrada en el régimen de sintonía
  const rendVol = (mot, n) =>
    Math.max(mot.etaVmax * (1 - 0.35 * Math.pow((n - mot.nEtaV) / mot.nEtaV, 2)), 0.55);

  // --------------------------------------------------------- integración RK4
  function integrar(geo, mot, cb, Qtot, m, p0, V0, n, th0, th1, opc) {
    const Sp = 2 * geo.S * n / 60;
    const soc = cb ? cb.soc : Infinity;
    const C1 = geo.Vd * mot.Tivc / (p0 * V0);
    const f = (t, p) => {
      const [V, dV] = volumen(geo, mot.rc, t);
      const dxb = cb ? quemado(t, cb)[1] : 0;
      const T = p * V / (m * AMB.R), g = gammaT(T);
      let dQw = 0;
      if (opc.woschni) {                                              // Woschni (1967)
        const pmot = p0 * Math.pow(V0 / V, mot.g);
        const C2 = t >= soc ? 3.24e-3 : 0;
        const w = Math.max(2.28 * Sp + C2 * C1 * (p - pmot), 0.1);
        const h = 3.26 * Math.pow(geo.B, -0.2) * Math.pow(p / 1e3, 0.8) * Math.pow(T, -0.55) * Math.pow(w, 0.8);
        const A = 2 * geo.Ap + 4 * V / geo.B;
        dQw = h * A * (T - opc.Tpared) / (6 * n);                    // [J/deg]
      }
      return (g - 1) / V * (Qtot * dxb - dQw) - g * p / V * dV;
    };
    const N = Math.round((th1 - th0) / opc.dth);
    const h = opc.dth / opc.sub;
    const th = new Float64Array(N + 1), p = new Float64Array(N + 1);
    th[0] = th0; p[0] = p0;
    let y = p0, t = th0;
    for (let i = 1; i <= N; i++) {
      for (let k = 0; k < opc.sub; k++) {
        const k1 = f(t, y), k2 = f(t + h / 2, y + h / 2 * k1);
        const k3 = f(t + h / 2, y + h / 2 * k2), k4 = f(t + h, y + h * k3);
        y += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        t += h;
      }
      th[i] = th0 + i * opc.dth; p[i] = y;
    }
    return { th, p };
  }

  // ------------------------------------------- autoignición (Livengood-Wu)
  // Retraso de iso-octano, He et al. (Combust. Flame 2005). p en atm, χO2 en %.
  // Validada en ~943-1027 K y 9-23 atm: fuera de ahí es extrapolación.
  const tauIgn = (p, T, phi, xO2) =>
    1.3e-4 * Math.pow(p / 101325, -1.05) * Math.pow(phi, -0.77) *
    Math.pow(100 * xO2, -1.41) * Math.exp(33700 / (1.987 * T));      // [ms]

  function socLW(mot, lam, p0, n, geo = GEO, opc = OPC) {
    const [V0] = volumen(geo, mot.rc, -180);
    const m = p0 * V0 / (AMB.R * mot.Tivc);
    const { th, p } = integrar(geo, mot, null, 0, m, p0, V0, n, -180, opc.finLW, opc);
    const dtms = 1000 / (6 * n), phi = 1 / lam;
    let I = 0, fPrev = null;
    for (let i = 0; i < th.length; i++) {
      const [V] = volumen(geo, mot.rc, th[i]);
      const fi = 1 / tauIgn(p[i], p[i] * V / (m * AMB.R), phi, mot.xO2);
      if (fPrev !== null) {
        const Inew = I + 0.5 * (fi + fPrev) * opc.dth * dtms;
        if (Inew >= 1) return th[i - 1] + (1 - I) / (Inew - I) * opc.dth;
        I = Inew;
      }
      fPrev = fi;
    }
    return null;                                                      // no autoinflama
  }

  const z50 = cb => Math.pow(Math.log(2) / cb.a, 1 / (cb.m + 1));

  // T_IVC que lleva el CA50 al objetivo (bisección: el SOC baja al subir T_IVC)
  function controlarCA50(mot, cb, lam, p0, n, geo = GEO, opc = OPC, Tlim = [330, 560], tol = 0.25) {
    const socObj = mot.ca50Obj - cb.dur * z50(cb);
    const soc = Ti => { const s = socLW(Object.assign({}, mot, { Tivc: Ti }), lam, p0, n, geo, opc); return s === null ? Infinity : s; };
    let [lo, hi] = Tlim;
    if (soc(hi) > socObj) return hi;
    if (soc(lo) < socObj) return lo;
    while (hi - lo > tol) {
      const mid = 0.5 * (lo + hi);
      if (soc(mid) > socObj) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  // ------------------------------------------------------------------- ciclo
  function ciclo(motIn, carga, n, opc = OPC, geo = GEO) {
    let mot = motIn;
    let { pIn, lam, cb, modo } = estrategia(mot, carga);
    const p0 = rendVol(mot, n) * pIn;
    if (mot.ctrl === 'hcci' && modo === 'Chispa') mot = Object.assign({}, mot, { Tivc: mot.TivcSI });
    if (cb.soc === null && mot.ca50Obj !== null && mot.ca50Obj !== undefined)
      mot = Object.assign({}, mot, { Tivc: controlarCA50(mot, cb, lam, p0, n, geo, opc) });

    const [V0] = volumen(geo, mot.rc, -180);
    const m = p0 * V0 / (AMB.R * mot.Tivc);
    const mf = m / (lam * mot.AFR);
    let Qtot = mot.etaComb * mf * mot.LHV, fallo = false;

    if (cb.soc === null) {
      const s = socLW(mot, lam, p0, n, geo, opc);
      if (s === null) { fallo = true; cb = null; Qtot = 0; modo = 'Fallo'; }
      else cb = Object.assign({}, cb, { soc: s });
    }

    const { th, p } = integrar(geo, mot, cb, Qtot, m, p0, V0, n, -180, 180, opc);
    const N = th.length;
    const V = new Float64Array(N), T = new Float64Array(N), Tb = new Array(N), hrr = new Float64Array(N);
    const lamFl = mot.llama === 'difusion' ? 1 : lam;                // llama de difusión: λ local ≈ 1
    const dTad = mot.etaComb * mot.LHV / ((1 + lamFl * mot.AFR) * 1450);
    const xO2b = Math.max(0.21 * (lamFl - 1) / lamFl, 0.01);
    let W = 0, nox = 0, pmax = 0, dpmax = -Infinity, Tmax = 0, ca50 = null, xPrev = 0, tasaPrev = 0;
    for (let i = 0; i < N; i++) {
      V[i] = volumen(geo, mot.rc, th[i])[0];
      T[i] = p[i] * V[i] / (m * AMB.R);
      const [xb, dxb] = quemado(th[i], cb);
      hrr[i] = Qtot * dxb;
      Tb[i] = xb >= 0.01 ? mot.Tivc * Math.pow(p[i] / p0, (mot.g - 1) / mot.g) + dTad : null;
      // Zeldovich con [O] en equilibrio (Heywood, ec. 11.8): d[NO]/dt ∝ T^-1/2·exp(-69090/T)·[O2]^1/2
      const tasa = Tb[i] ? (lamFl / lam) * xb * Math.pow(Tb[i], -0.5) * Math.exp(-69090 / Tb[i]) * Math.sqrt(xO2b) : 0;
      if (i > 0) {
        W += 0.5 * (p[i] + p[i - 1]) * (V[i] - V[i - 1]);
        nox += 0.5 * (tasa + tasaPrev) * opc.dth;
        dpmax = Math.max(dpmax, (p[i] - p[i - 1]) / opc.dth);
        if (ca50 === null && cb && xPrev < 0.5 && xb >= 0.5)
          ca50 = th[i - 1] + (0.5 - xPrev) / (xb - xPrev) * opc.dth;
      }
      pmax = Math.max(pmax, p[i]); Tmax = Math.max(Tmax, T[i]);
      xPrev = xb; tasaPrev = tasa;
    }
    const pmep = AMB.pEsc - pIn, imepG = W / geo.Vd;
    return {
      motor: motIn.nombre, modo, carga, n, lam, Tivc: mot.Tivc, fallo,
      soc: cb ? cb.soc : null, ca50,
      th, V, p, T, Tb, hrr,
      imepG, pmep, imepN: imepG - pmep, etaI: (W - pmep * geo.Vd) / (mf * mot.LHV),
      pmax, dpmax, Tmax, nox: nox / (6 * n)
    };
  }

  // Par [Nm] y potencia [CV] con FMEP de Chen-Flynn (coeficientes orientativos)
  function parEfectivo(r, geo = GEO) {
    const cm = 2 * geo.S * r.n / 60;
    const fmep = (0.25 + 0.005 * r.pmax / 1e5 + 0.06 * cm + 0.0009 * cm * cm) * 1e5;
    const par = (r.imepN - fmep) * geo.Vd * geo.nCil / (4 * Math.PI);
    return { par, pot: par * 2 * Math.PI * r.n / 60 / 735.5, fmep };
  }

  const api = { AMB, GEO, MOTORES, OPC, volumen, wiebe, quemado, estrategia, rendVol,
                tauIgn, socLW, controlarCA50, z50, ciclo, parEfectivo };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.Motor0D = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
