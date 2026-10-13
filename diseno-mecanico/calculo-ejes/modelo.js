/* Ultimate Shaft Solver v1.5 · motor de cálculo
 * Port en JavaScript de ultimate_shaft_solver.py (Shigley / Norton). Mismas
 * ecuaciones; brentq se sustituye por Brent propio y minimize_scalar por
 * sección áurea acotada.
 *
 * Unidades internas: mm, N, N·mm, MPa, kg/mm³.
 */
(function (root) {
  'use strict';

  const m = (nombre, fam, Sy, Sut, E, rho, Se_prime, N_e, ka, src) =>
    ({ nombre, family: fam, S_y: Sy, S_ut: Sut, E, rho, Se_prime, N_e, ka, src });

  const MATERIALS = {
    c45:      m('Acero C45 normalizado (+N)', 'steel', 340, 620, 210000, 7.85e-6, null, 1e6, 'machined', 'EN 10083-2, d ≤ 16 mm, mín.'),
    '42crmo4':m('Acero 42CrMo4 bonificado (+QT)', 'steel', 750, 1000, 210000, 7.85e-6, null, 1e6, 'machined', 'EN 10083-3, 16 < d ≤ 40 mm, mín.'),
    bolt_88:  m('Tornillo clase 8.8', 'bolt', 640, 800, 210000, 7.85e-6, null, 1e6, 'machined', 'ISO 898-1, nominal'),
    bolt_109: m('Tornillo clase 10.9', 'bolt', 940, 1040, 210000, 7.85e-6, null, 1e6, 'machined', 'ISO 898-1, nominal'),
    bolt_129: m('Tornillo clase 12.9', 'bolt', 1100, 1220, 210000, 7.85e-6, null, 1e6, 'machined', 'ISO 898-1, nominal'),
    aisi304:  m('Inoxidable AISI 304 recocido', 'stainless', 205, 515, 193000, 8.00e-6, 0.45 * 515, 1e7, 'machined', "ASTM A276, mín.; Se' ≈ 0,45·Sut"),
    al6061:   m('Aluminio 6061-T6', 'aluminum', 276, 310, 68900, 2.70e-6, 96.5, 5e8, 0.85, "ASM, típicos; Sf' R.R. Moore a 5e8"),
    al6082:   m('Aluminio 6082-T6', 'aluminum', 250, 290, 70000, 2.70e-6, 95.0, 5e8, 0.85, "EN 755-2, mín.; Sf' típico a 5e8"),
    al7020:   m('Aluminio 7020-T6', 'aluminum', 290, 350, 70000, 2.78e-6, null, 5e8, 0.85, "EN 755-2, mín.; Sf' regla de Norton"),
    al7075:   m('Aluminio 7075-T6', 'aluminum', 503, 572, 71700, 2.81e-6, 159.0, 5e8, 0.85, "ASM, típicos; Sf' R.R. Moore a 5e8"),
    ti64:     m('Titanio Ti-6Al-4V (grado 5) recocido', 'titanium', 880, 950, 113800, 4.43e-6, 510.0, 1e7, 'machined', "ASM, típicos; Se' probeta lisa R = −1 a 1e7"),
  };

  const FATIGUE_OPTS = { kd: 1.0, ke: 0.814, kx: 1.0 };   // ke 0,814 = 99 % de fiabilidad
  const KC_AXIAL = 0.85;
  const S1000_RATIO_AL = 0.9;

  // ------------------------------------------------------------ geometría
  function geometria(D, t) {
    const d_int = (t <= 0 || t >= D / 2) ? 0 : D - 2 * t;
    const area = Math.PI / 4 * (D ** 2 - d_int ** 2);
    const W = Math.PI * (D ** 4 - d_int ** 4) / (32 * D);
    const J = Math.PI * (D ** 4 - d_int ** 4) / 32;
    return { d_int, area, W, J };
  }

  function factor_kb(d) {                       // Shigley ec. 6-20 (mm), continuo y ≤ 1
    if (d <= 7.62) return 1.0;
    if (d <= 51) return 1.24 * d ** -0.107;
    return 1.51 * d ** -0.157;
  }

  function factor_f(Sut, family) {              // fracción de Sut a 1000 ciclos (Shigley fig. 6-18)
    if (family !== 'steel' && family !== 'bolt') return S1000_RATIO_AL;
    let k = Sut / 6.894757;
    if (k < 70) return 0.9;
    k = Math.min(k, 200);
    return 1.06 - 2.8e-3 * k + 6.9e-6 * k * k;
  }

  function limite_fatiga(mat, D, gira, solo_axial) {
    let S_prime;
    if (mat.Se_prime !== null) S_prime = mat.Se_prime;
    else if (['steel', 'bolt', 'stainless'].includes(mat.family)) S_prime = mat.S_ut <= 1400 ? 0.5 * mat.S_ut : 700;
    else if (mat.family === 'aluminum') S_prime = mat.S_ut < 330 ? 0.4 * mat.S_ut : 130;   // Norton, a 5e8
    else S_prime = 0.5 * mat.S_ut;

    let ka = mat.ka === 'machined' ? 4.51 * mat.S_ut ** -0.265 : mat.ka;
    ka = Math.min(ka, 1.0);
    const kb = solo_axial ? 1.0 : factor_kb(gira ? D : 0.37 * D);
    const kc = solo_axial ? KC_AXIAL : 1.0;     // combinadas: el 0,85 axial va dentro del VM
    const { kd, ke, kx } = FATIGUE_OPTS;
    const Se = S_prime * ka * kb * kc * kd * ke * kx;
    return { S_prime, factors: { ka, kb, kc, kd, ke, kx }, Se };
  }

  // cargas = { flex: [pico, media, alt] (N·mm), axil: [..] (N), tor: [..] (N·mm) }
  function calcular(D, t, cargas, kf_f, kf_a, kf_t, matKey, gira = true) {
    const mat = MATERIALS[matKey];
    const { d_int, area, W, J } = geometria(D, t);
    const r = D / 2;
    const [, M_m, M_a] = cargas.flex, [, F_m, F_a] = cargas.axil, [, T_m, T_a] = cargas.tor;

    // Flexión: la fibra crítica es la traccionada -> |M_m|. Axial: el signo se conserva.
    const sb_m = Math.abs(M_m) / W, sb_a = M_a / W;
    const sx_ax_m = F_m / area, sx_ax_a = F_a / area;
    const t_m = T_m * r / J, t_a = T_a * r / J;

    const solo_axial = M_a === 0 && M_m === 0 && T_a === 0 && T_m === 0 && (F_a !== 0 || F_m !== 0);

    // Fatiga (con Kf). Axial alternante dividido entre 0,85 (Shigley 7-7)
    const kcVM = solo_axial ? 1.0 : KC_AXIAL;
    const sx_a = kf_f * Math.abs(sb_a) + kf_a * Math.abs(sx_ax_a) / kcVM;
    const sx_m = kf_f * sb_m + kf_a * sx_ax_m;
    const vm_a = Math.sqrt(sx_a ** 2 + 3 * (kf_t * t_a) ** 2);
    const compresion = sx_m < 0;
    const vm_m = Math.sqrt((compresion ? 0 : sx_m) ** 2 + 3 * (kf_t * t_m) ** 2);

    // Estática (sin Kf, material dúctil): pico nominal
    const sx_peak = Math.abs(sb_m) + Math.abs(sb_a) + Math.abs(sx_ax_m) + Math.abs(sx_ax_a);
    const t_peak = Math.abs(t_m) + Math.abs(t_a);
    const vm_peak = Math.sqrt(sx_peak ** 2 + 3 * t_peak ** 2);

    const { S_prime, factors, Se } = limite_fatiga(mat, D, gira, solo_axial);
    const Sy = mat.S_y, Sut = mat.S_ut;

    const fos_est = vm_peak > 0 ? Sy / vm_peak : Infinity;
    const vm_m_full = Math.sqrt(sx_m ** 2 + 3 * (kf_t * t_m) ** 2);
    const fos_lang = (vm_a + vm_m_full) > 0 ? Sy / (vm_a + vm_m_full) : Infinity;
    const tg = vm_a / Se + vm_m / Sut;
    const fos_good = tg > 0 ? 1 / tg : Infinity;
    const ta = (vm_a / Se) ** 2 + (vm_m / Sy) ** 2;
    const fos_asme = ta > 0 ? 1 / Math.sqrt(ta) : Infinity;

    // Vida (Basquin con Goodman)
    const N_e = mat.N_e;
    const life = { type: null, N: null };
    const S_1000 = factor_f(Sut, mat.family) * Sut;
    let sig_rev;
    if (vm_m >= Sut) { sig_rev = NaN; life.type = 'static'; }
    else {
      sig_rev = vm_a / (1 - vm_m / Sut);
      if (sig_rev <= Se) life.type = 'inf';
      else if (sig_rev >= S_1000) life.type = 'low';
      else {
        const b = Math.log10(S_1000 / Se) / Math.log10(1e3 / N_e);
        const a = S_1000 / 1e3 ** b;
        life.type = 'finite'; life.N = (sig_rev / a) ** (1 / b);
      }
    }

    return {
      d_int, area, vm_peak, vm_m, vm_a, compresion, Se, S_prime, N_e, factors,
      fos: { static: fos_est, langer: fos_lang, goodman: fos_good, asme: fos_asme },
      sig_rev, life, sx_peak, t_peak, S_1000, S_y: Sy, S_ut: Sut, family: mat.family, vm_m_full,
      masa_por_mm: area * mat.rho * 1000,          // g/mm
    };
  }

  function fos_criterio(D, t, a) {
    const r = calcular(D, t, a.cargas, a.kf_f, a.kf_a, a.kf_t, a.matKey, a.gira).fos;
    return a.modo === 'stat' ? r.static : Math.min(r.goodman, r.asme);
  }

  // Brent (equivalente a scipy.optimize.brentq)
  function brent(f, a, b, xtol = 1e-4, maxit = 200) {
    let fa = f(a), fb = f(b);
    if (fa * fb > 0) return null;
    let c = a, fc = fa, d = b - a, e = d;
    for (let i = 0; i < maxit; i++) {
      if (fb * fc > 0) { c = a; fc = fa; d = e = b - a; }
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      const tol = 2 * Number.EPSILON * Math.abs(b) + xtol / 2, mm = (c - b) / 2;
      if (Math.abs(mm) <= tol || fb === 0) return b;
      if (Math.abs(e) >= tol && Math.abs(fa) > Math.abs(fb)) {
        let p, q, s = fb / fa;
        if (a === c) { p = 2 * mm * s; q = 1 - s; }
        else { const qq = fa / fc, rr = fb / fc; p = s * (2 * mm * qq * (qq - rr) - (b - a) * (rr - 1)); q = (qq - 1) * (rr - 1) * (s - 1); }
        if (p > 0) q = -q; else p = -p;
        if (2 * p < Math.min(3 * mm * q - Math.abs(tol * q), Math.abs(e * q))) { e = d; d = p / q; }
        else { d = mm; e = d; }
      } else { d = mm; e = d; }
      a = b; fa = fb;
      b += Math.abs(d) > tol ? d : (mm > 0 ? tol : -tol);
      fb = f(b);
    }
    return b;
  }

  function resolverD(espDeD, Dmin, Dmax, args, target) {
    const g = D => fos_criterio(D, espDeD(D), args) - target;
    if (g(Dmin) >= 0) return Dmin;
    if (g(Dmax) < 0) return null;
    return brent(g, Dmin, Dmax);
  }

  function aurea(f, lo, hi, xtol = 1e-4) {
    const r = (Math.sqrt(5) - 1) / 2;
    let a = lo, b = hi, c = b - r * (b - a), d = a + r * (b - a), fc = f(c), fd = f(d);
    while (b - a > xtol) {
      if (fc < fd) { b = d; d = c; fd = fc; c = b - r * (b - a); fc = f(c); }
      else { a = c; c = d; fc = fd; d = a + r * (b - a); fd = f(d); }
    }
    const x = (a + b) / 2;
    return { x, fun: f(x) };
  }

  /* geoModo: 'fixed_t' | 'ratio' | 'mass'. Devuelve {D, t} o null. */
  function optimizar(geoModo, tIni, DIni, args, target, { Dmax = 1000, tmin = 1, DTmax = 20 } = {}) {
    const macizo = tIni <= 0 || tIni >= DIni / 2;
    if (geoModo === 'fixed_t' || (geoModo === 'ratio' && macizo)) {
      if (macizo) { const D = resolverD(D => D / 2, 0.5, Dmax, args, target); return D === null ? null : { D, t: D / 2 }; }
      const D = resolverD(() => tIni, 2 * tIni + 0.5, Dmax, args, target);
      return D === null ? null : { D, t: tIni };
    }
    if (geoModo === 'ratio') {
      const k = (DIni - 2 * tIni) / DIni;
      const esp = D => D * (1 - k) / 2;
      const D = resolverD(esp, 0.5, Dmax, args, target);
      return D === null ? null : { D, t: esp(D) };
    }
    // 'mass': barrido logarítmico en t; para cada t el D mínimo; menor área
    const area = (D, t) => geometria(D, t).area;
    const cand = [];
    const Ds = resolverD(D => D / 2, 0.5, Dmax, args, target);
    if (Ds !== null) cand.push([area(Ds, Ds / 2), Ds, Ds / 2]);
    const thi = Dmax / 2, n = 120;
    for (let i = 0; i <= n; i++) {
      const t = tmin * (thi / tmin) ** (i / n);
      const Dlo = Math.max(2 * t + 0.5, 0.5);
      if (Dlo > Dmax) break;
      const Dhi = Math.min(Dmax, DTmax * t);
      if (Dlo > Dhi) continue;
      const D = resolverD(() => t, Dlo, Dhi, args, target);
      if (D !== null && D / t <= DTmax + 1e-9) cand.push([area(D, t), D, t]);
    }
    if (!cand.length) return null;
    cand.sort((p, q) => p[0] - q[0]);
    let [Abest, D, t] = cand[0];
    if (2 * t < D) {
      const paso = (thi / tmin) ** (1 / n);
      const areaDeT = tt => {
        const Dlo = 2 * tt + 0.5, Dhi = Math.min(Dmax, DTmax * tt);
        if (Dlo > Dhi) return Infinity;
        const Dt = resolverD(() => tt, Dlo, Dhi, args, target);
        return Dt === null ? Infinity : area(Dt, tt);
      };
      const r = aurea(areaDeT, Math.max(tmin, t / paso), t * paso);
      if (r.fun < Abest) {
        t = r.x;
        D = resolverD(() => t, 2 * t + 0.5, Math.min(Dmax, DTmax * t), args, target);
      }
    }
    return { D, t };
  }

  const api = { MATERIALS, FATIGUE_OPTS, geometria, factor_kb, factor_f, limite_fatiga, calcular, fos_criterio, optimizar };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.USS = api;
})(this);
