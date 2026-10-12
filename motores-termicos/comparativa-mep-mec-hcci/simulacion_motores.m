% =========================================================================
% Modelo 0D de ciclo cerrado: comparativa MEP / MEC / HCCI  (versión MATLAB)
% -------------------------------------------------------------------------
% Equivalente a la versión Python del repositorio (motor0d/).
% - Monozona, gas ideal, gamma(T) en la ecuación de energía
% - Liberación de calor: Wiebe (MEP/HCCI), Wiebe doble de Watson (MEC)
% - Pérdidas a pared: Woschni (1967)
% - Control de carga: cantidad (mariposa + bombeo) vs calidad (lambda)
% - HCCI: autoignición por Livengood-Wu (He et al., 2005) y control de CA50
%   en lazo cerrado actuando sobre la T en cierre de admisión (T_IVC)
% - NOx: índice cinético tipo Zeldovich sobre T estimada de quemados
% - Par/potencia: IMEP neta - FMEP (Chen-Flynn)
% Requiere MATLAB R2020a+ (funciones locales en script, tiledlayout, exportgraphics).
% Ejecutar desde la carpeta matlab/. Validado frente a la versión Python (Octave 8.4, modelo).
% =========================================================================
clear; clc; close all;

%% 1. Configuración
amb  = struct('p', 1.013e5, 'T', 298, 'R', 287, 'pExh', 1.05*1.013e5);
geo  = defGeometria();   % misma geometría para los tres: comparativa justa
mot  = defMotores();
opts = struct('woschni', true, 'Tw', 450, 'dth', 0.5, 'finLW', 30, 'exportar', false);

cargaRef = 0.40;          % dentro de la ventana HCCI
nRef     = 2500;          % [rpm]
cargas   = 0.15:0.05:1.0;
rpm      = 1000:250:6500;
nM       = numel(mot);

%% 2. Ciclo de referencia
for k = 1:nM
    ref(k) = cicloCerrado(mot(k), geo, amb, opts, cargaRef, nRef); %#ok<SAGROW>
end

fprintf('\nPunto de referencia: carga %.0f %%, %d rpm\n', cargaRef*100, nRef);
fprintf('%-6s %-6s %6s %7s %7s %7s %6s %7s %10s %7s\n', 'Motor', 'Modo', 'lambda', ...
        'T_IVC', 'CA50', 'IMEPn', 'eta_i', 'pmax', 'dp/dth', 'Tmax');
for k = 1:nM
    r = ref(k);
    fprintf('%-6s %-6s %6.2f %7.0f %+7.1f %7.2f %6.1f %7.1f %10.2f %7.0f\n', ...
        mot(k).nombre, r.modo, r.lambda, r.Tivc, r.ca50, r.imepN/1e5, r.etaI*100, ...
        r.pmax/1e5, r.dpmax/1e5, r.Tmax);
end

%% 3. Barrido de carga
[L, IM, ET, NX] = deal(nan(nM, numel(cargas)));
for k = 1:nM
    for j = 1:numel(cargas)
        r = cicloCerrado(mot(k), geo, amb, opts, cargas(j), nRef);
        L(k,j) = r.lambda;  IM(k,j) = r.imepN;
        ET(k,j) = r.etaI;   NX(k,j) = r.noxIdx;
    end
end
NX = NX / NX(1, end);     % normalizado: MEP a plena carga = 1

%% 4. Curvas a plena carga
[Tq, Pw] = deal(nan(nM, numel(rpm)));
for k = 1:nM
    for j = 1:numel(rpm)
        if rpm(j) > mot(k).nMax, continue; end
        r = cicloCerrado(mot(k), geo, amb, opts, 1.0, rpm(j));
        [Tq(k,j), Pw(k,j)] = parEfectivo(r, geo, rpm(j));
    end
end

%% 5. Control HCCI: lazo abierto (CA50 vs T_IVC) y lazo cerrado (T_IVC vs carga)
h      = mot(3);
rpmsH  = [1500 2500 4000];
TivcV  = 400:2.5:520;
cargasH = cargas(cargas <= 1/h.lambdaMin + 1e-9);
cb     = h.comb;
z50    = (log(2)/cb.a)^(1/(cb.m + 1));
[pInH, lamH] = estrategiaCarga(h, amb, cargaRef);
CA50la = nan(numel(rpmsH), numel(TivcV));
Treq   = nan(numel(rpmsH), numel(cargasH));
for i = 1:numel(rpmsH)
    p0 = rendVol(h, rpmsH(i))*pInH;
    for j = 1:numel(TivcV)
        hj = h;  hj.Tivc = TivcV(j);
        s = socLW(hj, lamH, p0, rpmsH(i), geo, amb, opts);
        if ~isempty(s), CA50la(i,j) = s + cb.dur*z50; end
    end
    for j = 1:numel(cargasH)
        r = cicloCerrado(h, geo, amb, opts, cargasH(j), rpmsH(i));
        Treq(i,j) = r.Tivc;
    end
end

%% 6. Figuras
f1 = figure('Name', 'Ciclo de referencia', 'Position', [80 80 1100 750]);
t1 = tiledlayout(f1, 2, 2, 'TileSpacing', 'compact');
title(t1, sprintf('Ciclo de referencia: carga %.0f %%, %d rpm', cargaRef*100, nRef));

nexttile; hold on; grid on;
for k = 1:nM
    plot(ref(k).V*1e6, ref(k).p/1e5, 'Color', mot(k).color, 'LineWidth', 1.8);
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('V [cm^3]'); ylabel('p [bar]'); title('Diagrama p-V (log-log)');
legend({mot.nombre}, 'Location', 'northeast');

nexttile; hold on; grid on;
for k = 1:nM
    plot(ref(k).th, ref(k).p/1e5, 'Color', mot(k).color, 'LineWidth', 1.8);
end
xlim([-60 90]); xlabel('\theta [deg]'); ylabel('p [bar]'); title('Presión en cámara');

nexttile; hold on; grid on;
for k = 1:nM
    plot(ref(k).th, ref(k).hrr, 'Color', mot(k).color, 'LineWidth', 1.8);
end
xlim([-40 90]); xlabel('\theta [deg]'); ylabel('dQ/d\theta [J/deg]');
title('Tasa de liberación de calor');

nexttile; hold on; grid on;
for k = 1:nM
    plot(ref(k).th, ref(k).T,  '-',  'Color', mot(k).color, 'LineWidth', 1.8);
    plot(ref(k).th, ref(k).Tb, '--', 'Color', mot(k).color, 'LineWidth', 1.2, ...
         'HandleVisibility', 'off');
end
yline(1800, ':k', 'NO térmico relevante (~1800 K)', 'HandleVisibility', 'off', ...
      'LabelHorizontalAlignment', 'left');
xlim([-40 90]); xlabel('\theta [deg]'); ylabel('T [K]');
title('T media (—) y de quemados estimada (- -)');

f2 = figure('Name', 'Barrido de carga', 'Position', [120 120 1100 750]);
t2 = tiledlayout(f2, 2, 2, 'TileSpacing', 'compact');
title(t2, sprintf('Barrido de carga a %d rpm', nRef));
datos  = {L, IM/1e5, ET*100, NX};
etiqY  = {'\lambda [-]', 'IMEP neta [bar]', '\eta_{ind} [%]', 'Índice NOx [-] (MEP PC = 1)'};
esLog  = [false false false true];
for i = 1:4
    nexttile; hold on; grid on;
    for k = 1:nM
        plot(cargas*100, datos{i}(k,:), 'Color', mot(k).color, 'LineWidth', 1.8);
    end
    if esLog(i), set(gca, 'YScale', 'log'); end
    xline(100/h.lambdaMin, '--k', 'Fin ventana HCCI', 'HandleVisibility', 'off');
    xlabel('Carga [%]'); ylabel(etiqY{i});
    if i == 1, legend({mot.nombre}, 'Location', 'northeast'); end
end

f3 = figure('Name', 'Plena carga', 'Position', [160 160 1100 420]);
t3 = tiledlayout(f3, 1, 2, 'TileSpacing', 'compact');
title(t3, 'Plena carga (HCCI en modo chispa; picado no modelado)');
nexttile; hold on; grid on;
for k = 1:nM, plot(rpm, Tq(k,:), 'Color', mot(k).color, 'LineWidth', 1.8); end
xlabel('n [rpm]'); ylabel('Par [Nm]');
nexttile; hold on; grid on;
for k = 1:nM, plot(rpm, Pw(k,:), 'Color', mot(k).color, 'LineWidth', 1.8); end
xlabel('n [rpm]'); ylabel('Potencia [CV]'); legend({mot.nombre}, 'Location', 'northwest');

f4 = figure('Name', 'Control HCCI', 'Position', [200 200 1100 420]);
t4 = tiledlayout(f4, 1, 2, 'TileSpacing', 'compact');
title(t4, 'Control de la autoignición HCCI (Livengood-Wu)');
secuencial = [0.62 0.85 0.75; 0.11 0.69 0.48; 0.04 0.42 0.29];   % un solo tono
etiq = arrayfun(@(n) sprintf('%d rpm', n), rpmsH, 'UniformOutput', false);
nexttile; hold on; grid on;
for i = 1:numel(rpmsH)
    plot(TivcV, CA50la(i,:), 'Color', secuencial(i,:), 'LineWidth', 1.8);
end
yline(h.ca50obj, ':k', sprintf('objetivo CA50 = %+.0f°', h.ca50obj), 'HandleVisibility', 'off');
xlabel('T en IVC [K]'); ylabel('CA50 [° tras PMS]'); legend(etiq);
title(sprintf('Lazo abierto (carga %.0f %%); sin curva = fallo', cargaRef*100));
nexttile; hold on; grid on;
for i = 1:numel(rpmsH)
    plot(cargasH*100, Treq(i,:), '-o', 'Color', secuencial(i,:), 'LineWidth', 1.8, ...
         'MarkerSize', 4, 'MarkerFaceColor', secuencial(i,:));
end
xlabel('Carga [%]'); ylabel('T en IVC [K]'); legend(etiq);
title(sprintf('Lazo cerrado: T_{IVC} para CA50 = %+.0f°', h.ca50obj));

if opts.exportar
    if ~exist('figs', 'dir'), mkdir('figs'); end
    exportgraphics(f1, fullfile('figs', 'ciclo_referencia.png'), 'Resolution', 200);
    exportgraphics(f2, fullfile('figs', 'barrido_carga.png'),    'Resolution', 200);
    exportgraphics(f3, fullfile('figs', 'plena_carga.png'),      'Resolution', 200);
    exportgraphics(f4, fullfile('figs', 'control_hcci.png'),     'Resolution', 200);
end

% =========================================================================
%  DEFINICIONES
% =========================================================================
function geo = defGeometria()
    geo.B = 0.086;  geo.S = 0.086;  geo.L = 0.145;  geo.nCil = 4;  % 2.0 L
    geo.Ap = pi*geo.B^2/4;
    geo.Vd = geo.Ap*geo.S;                                          % por cilindro
end

function mot = defMotores()
    % Campos comunes; los específicos de HCCI quedan vacíos en MEP y MEC.
    % g: gamma de compresión de la carga fresca (arrastre y T de no quemados)

    % --- MEP: control por cantidad (mariposa), lambda = 1, premezcla
    mot(1).nombre = 'MEP';  mot(1).color = [0.165 0.471 0.839];   % #2a78d6
    mot(1).rc = 10.5;  mot(1).g = 1.32;  mot(1).Tivc = 320;
    mot(1).LHV = 43.0e6;  mot(1).AFRst = 14.7;  mot(1).etaComb = 0.97;
    mot(1).ctrl = 'cantidad';  mot(1).llama = 'premezcla';  mot(1).lambdaPC = 1.0;
    mot(1).comb = struct('tipo', 'wiebe', 'soc', -22, 'dur', 50, 'a', 5, 'm', 2);
    mot(1).nMax = 6500;  mot(1).nEtaV = 4500;  mot(1).etaVmax = 0.92;
    mot(1).xO2 = 0.21;

    % --- MEC: control por calidad, fase premezclada + difusión
    mot(2).nombre = 'MEC';  mot(2).color = [0.922 0.408 0.204];   % #eb6834
    mot(2).rc = 17.5;  mot(2).g = 1.35;  mot(2).Tivc = 330;
    mot(2).LHV = 42.6e6;  mot(2).AFRst = 14.5;  mot(2).etaComb = 0.98;
    mot(2).ctrl = 'calidad';  mot(2).llama = 'difusion';  mot(2).lambdaPC = 1.35; % humos
    mot(2).comb = struct('tipo', 'doble', 'soc', -6, 'dur', 8, 'a', 5, 'm', 2, ...
                         'beta', 0.15, 'dur2', 55, 'm2', 0.9);
    mot(2).nMax = 4500;  mot(2).nEtaV = 2200;  mot(2).etaVmax = 0.90;
    mot(2).xO2 = 0.21;

    % --- HCCI: pobre sin mariposa hasta lambdaMin; por encima, modo chispa
    mot(3).nombre = 'HCCI';  mot(3).color = [0.106 0.686 0.478];  % #1baf7a
    mot(3).rc = 15;  mot(3).g = 1.36;
    mot(3).Tivc = 440;       % valor inicial; con ca50obj lo fija el control
    mot(3).LHV = 43.0e6;  mot(3).AFRst = 14.7;  mot(3).etaComb = 0.95;  % más HC/CO
    mot(3).ctrl = 'hcci';  mot(3).llama = 'premezcla';  mot(3).lambdaPC = 1.0;
    mot(3).lambdaMin = 2.0;  % límite de ringing (dp/dtheta)
    % soc = NaN: el inicio de combustión lo calcula Livengood-Wu
    mot(3).comb   = struct('tipo', 'wiebe', 'soc', NaN, 'dur', 12, 'a', 5, 'm', 2);
    mot(3).combSI = struct('tipo', 'wiebe', 'soc', -8,  'dur', 55, 'a', 5, 'm', 2);
    mot(3).TivcSI = 330;     % modo chispa: sin residuales calientes
    mot(3).nMax = 6500;  mot(3).nEtaV = 4000;  mot(3).etaVmax = 0.88;
    mot(3).xO2 = 0.19;       % dilución por residuales
    mot(3).ca50obj = 6;      % control de fase: CA50 ~ 6 deg tras PMS ([] = lazo abierto)
end

% =========================================================================
%  MODELO
% =========================================================================
function [pIn, lambda, cb, modo] = estrategiaCarga(mot, amb, carga)
    % carga = masa de combustible / masa máxima del motor (lambdaPC, sin mariposa)
    cb = mot.comb;  modo = mot.nombre;
    switch mot.ctrl
        case 'cantidad'                         % mariposa, lambda fija
            lambda = 1;  pIn = carga*amb.p;
        case 'calidad'                          % sin mariposa, lambda variable
            lambda = mot.lambdaPC/carga;  pIn = amb.p;
            cb.beta = 0.10 + 0.35*(1 - carga);  % más premezcla a baja carga
            cb.dur2 = 30 + 30*carga;            % difusión más larga a alta carga
        case 'hcci'
            if carga <= 1/mot.lambdaMin         % ventana HCCI (mezcla pobre)
                lambda = 1/carga;  pIn = amb.p;  modo = 'HCCI';
            else                                % transición a chispa, lambda = 1
                lambda = 1;  pIn = carga*amb.p;  cb = mot.combSI;  modo = 'SI';
            end
    end
end

function etaV = rendVol(mot, n)
    % Curva SUPUESTA (no calculada): parábola centrada en el régimen de sintonía
    etaV = max(mot.etaVmax*(1 - 0.35*((n - mot.nEtaV)/mot.nEtaV).^2), 0.55);
end

function r = cicloCerrado(mot, geo, amb, opts, carga, n)
    [pIn, lambda, cb, modo] = estrategiaCarga(mot, amb, carga);
    th = (-180:opts.dth:180)';              % IVC = EVO = PMI (ciclo cerrado ideal)
    V  = volumen(geo, mot.rc, th);
    p0 = rendVol(mot, n)*pIn;               % presión en cierre de admisión
    if strcmp(modo, 'SI'), mot.Tivc = mot.TivcSI; end

    % --- Autoignición: control de CA50 sobre T_IVC y Livengood-Wu
    fallo = false;
    if isnan(cb.soc)
        if ~isempty(mot.ca50obj)
            mot.Tivc = controlarCA50(mot, cb, lambda, p0, n, geo, amb, opts);
        end
        soc = socLW(mot, lambda, p0, n, geo, amb, opts);
        if isempty(soc)
            fallo = true;  cb.soc = Inf;    % no autoinflama: ciclo arrastrado
            modo  = 'fallo';
        else
            cb.soc = soc;
        end
    end

    m  = p0*V(1)/(amb.R*mot.Tivc);          % masa atrapada
    mf = m/(lambda*mot.AFRst);
    Qtot = mot.etaComb*mf*mot.LHV*(~fallo);

    p = integrar(geo, mot, cb, Qtot, m, p0, V(1), n, th, opts, amb);
    T = p.*V/(m*amb.R);
    [xb, dxb] = quemado(th, cb);
    W = trapz(V, p);                        % trabajo bruto de ciclo cerrado [J]

    r.th = th;  r.V = V;  r.p = p;  r.T = T;  r.hrr = Qtot*dxb;
    r.lambda = lambda;  r.modo = modo;  r.Tivc = mot.Tivc;  r.soc = cb.soc;
    r.ca50 = NaN;
    if ~fallo, r.ca50 = interp1(xb(xb > 1e-6 & xb < 1 - 1e-6), ...
                                th(xb > 1e-6 & xb < 1 - 1e-6), 0.5); end
    r.pmax = max(p);  r.Tmax = max(T);  r.dpmax = max(diff(p))/opts.dth;
    r.imepG = W/geo.Vd;
    r.pmep  = amb.pExh - pIn;               % bombeo simplificado
    r.imepN = r.imepG - r.pmep;
    r.etaI  = (W - r.pmep*geo.Vd)/(mf*mot.LHV);

    % --- T de quemados estimada y NOx (indicador cualitativo, no ppm)
    % En difusión la llama arde a lambda local ~1 aunque la global sea pobre
    if strcmp(mot.llama, 'difusion'), lamFl = 1; else, lamFl = lambda; end
    Tu   = mot.Tivc*(p/p0).^((mot.g - 1)/mot.g);          % no quemados, isentrópica
    dTad = mot.etaComb*mot.LHV/((1 + lamFl*mot.AFRst)*1450);
    Tb   = Tu + dTad;  Tb(xb < 0.01) = NaN;               % sin disociación: sobreestima
    r.Tb = Tb;
    xO2b = max(0.21*(lamFl - 1)/lamFl, 0.01);             % suelo ~ O2 por disociación
    tasa = (lamFl/lambda)*xb.*exp(-38370./Tb)*sqrt(xO2b);  tasa(isnan(tasa)) = 0;
    r.noxIdx = trapz(th, tasa)/(6*n);                     % integral temporal
end

function p = integrar(geo, mot, cb, Qtot, m, p0, V0, n, th, opts, amb)
    c = struct('geo', geo, 'rc', mot.rc, 'g', mot.g, 'R', amb.R, 'm', m, ...
               'comb', cb, 'Qtot', Qtot, 'woschni', opts.woschni, 'Tw', opts.Tw, ...
               'n', n, 'p0', p0, 'V0', V0, 'T0', mot.Tivc, 'Sp', 2*geo.S*n/60);
    odeOpt = odeset('MaxStep', 0.5, 'RelTol', 1e-6, 'AbsTol', 1);
    [~, p] = ode45(@(t, y) dpdth(t, y, c), th, p0, odeOpt);
end

function dp = dpdth(th, p, c)
    [V, dV] = volumen(c.geo, c.rc, th);
    [~, dxb] = quemado(th, c.comb);
    T  = p*V/(c.m*c.R);
    gT = 1.392 - 8.13e-5*T;          % gamma(T): con gamma cte. la T de quemados se dispara
    dQw = 0;
    if c.woschni
        pmot = c.p0*(c.V0/V)^c.g;
        C2   = 3.24e-3*(th >= c.comb.soc);
        w    = max(2.28*c.Sp + C2*c.geo.Vd*c.T0/(c.p0*c.V0)*(p - pmot), 0.1);
        h    = 3.26*c.geo.B^-0.2*(p/1e3)^0.8*T^-0.55*w^0.8;
        A    = 2*c.geo.Ap + 4*V/c.geo.B;
        dQw  = h*A*(T - c.Tw)/(6*c.n);                    % [J/deg]
    end
    dp = (gT - 1)/V*(c.Qtot*dxb - dQw) - gT*p/V*dV;
end

% -------------------------------------------------------------------------
%  Autoignición HCCI
% -------------------------------------------------------------------------
function tau = tauIgn(p, T, phi, xO2)
    % Retraso de autoignición [ms], iso-octano (He et al., Combust. Flame 2005)
    % p en atm, xO2 en % molar, Ru = 1.987 cal/mol K. Validada ~943-1027 K,
    % 9-23 atm, phi 0.25-1: fuera de ese rango es extrapolación.
    tau = 1.3e-4*(p/101325).^-1.05*phi^-0.77*(100*xO2)^-1.41.*exp(33700./(1.987*T));
end

function soc = socLW(mot, lambda, p0, n, geo, amb, opts)
    % Inicio de combustión [deg] por Livengood-Wu, o [] si no autoinflama.
    % Hasta el encendido no hay liberación de calor: la traza arrastrada es la real.
    th = (-180:opts.dth:opts.finLW)';
    V  = volumen(geo, mot.rc, th);
    m  = p0*V(1)/(amb.R*mot.Tivc);
    cb = mot.comb;  cb.soc = Inf;           % sin combustión
    p  = integrar(geo, mot, cb, 0, m, p0, V(1), n, th, opts, amb);
    T  = p.*V/(m*amb.R);
    I  = cumtrapz(th, 1./tauIgn(p, T, 1/lambda, mot.xO2))*1000/(6*n);  % ms/deg
    k  = find(I >= 1, 1);
    if isempty(k)
        soc = [];
    else
        soc = interp1(I(k-1:k), th(k-1:k), 1);
    end
end

function Tivc = controlarCA50(mot, cb, lambda, p0, n, geo, amb, opts)
    % T_IVC [K] que lleva el CA50 al objetivo (bisección; el SOC baja al subir T_IVC).
    % Representa el control real del HCCI: residuales internos (NVO) o EGR caliente.
    z50    = (log(2)/cb.a)^(1/(cb.m + 1));
    socObj = mot.ca50obj - cb.dur*z50;
    socDe  = @(Ti) socConTivc(mot, Ti, lambda, p0, n, geo, amb, opts);
    lo = 330;  hi = 560;  tol = 0.25;
    if socDe(hi) > socObj, Tivc = hi; return; end
    if socDe(lo) < socObj, Tivc = lo; return; end
    while hi - lo > tol
        mid = 0.5*(lo + hi);
        if socDe(mid) > socObj, lo = mid; else, hi = mid; end
    end
    Tivc = 0.5*(lo + hi);
end

function s = socConTivc(mot, Ti, lambda, p0, n, geo, amb, opts)
    mot.Tivc = Ti;
    s = socLW(mot, lambda, p0, n, geo, amb, opts);
    if isempty(s), s = Inf; end                 % fallo de encendido
end

% -------------------------------------------------------------------------
%  Utilidades
% -------------------------------------------------------------------------
function [V, dV] = volumen(geo, rc, th)
    % Biela-manivela exacto, theta = 0 en PMS. dV en [m^3/deg]
    t = deg2rad(th);  a = geo.S/2;  R = geo.L/a;
    Vc   = geo.Vd/(rc - 1);
    raiz = sqrt(R^2 - sin(t).^2);
    V  = Vc + geo.Ap*a*(R + 1 - cos(t) - raiz);
    dV = geo.Ap*a*sin(t).*(1 + cos(t)./raiz)*pi/180;
end

function [xb, dxb] = quemado(th, cb)
    switch cb.tipo
        case 'wiebe'
            [xb, dxb] = wiebe(th, cb.soc, cb.dur, cb.a, cb.m);
        case 'doble'   % Watson: premezcla (beta) + difusión
            [x1, d1] = wiebe(th, cb.soc, cb.dur,  cb.a, cb.m);
            [x2, d2] = wiebe(th, cb.soc, cb.dur2, cb.a, cb.m2);
            xb  = cb.beta*x1 + (1 - cb.beta)*x2;
            dxb = cb.beta*d1 + (1 - cb.beta)*d2;
    end
end

function [xb, dxb] = wiebe(th, soc, dur, a, m)
    z   = max((th - soc)/dur, 0);
    xb  = 1 - exp(-a*z.^(m + 1));
    dxb = a*(m + 1)/dur*z.^m.*exp(-a*z.^(m + 1));         % [1/deg]
end

function [Tq, P_CV] = parEfectivo(r, geo, n)
    cm   = 2*geo.S*n/60;
    % Chen-Flynn, coeficientes ORIENTATIVOS (calibrar con datos reales)
    fmep = (0.25 + 0.005*r.pmax/1e5 + 0.06*cm + 0.0009*cm^2)*1e5;
    bmep = r.imepN - fmep;
    Tq   = bmep*geo.Vd*geo.nCil/(4*pi);                   % 4T
    P_CV = Tq*2*pi*n/60/735.5;
end
