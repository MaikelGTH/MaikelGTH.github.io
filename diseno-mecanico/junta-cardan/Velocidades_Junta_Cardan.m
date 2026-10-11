%% CINEMÁTICA DE DOBLE JUNTA CARDAN (Z, W Y ESPACIAL)
clear; clc; close all;

%% PARÁMETROS
omega_1    = 100;   % Velocidad angular del motor (rad/s)
delta1_deg = 15;    % Ángulo de flexión de la 1ª junta (grados)
delta2_deg = 15;    % Ángulo de flexión de la 2ª junta (semiángulo del cono)

% Ángulo azimutal (beta): posición del Eje 3 en la base del cono
%   180° = Configuración en Z (plano XY, ejes 1 y 3 paralelos)
%   0°   = Configuración en W (plano XY)
%   90°  = Flexión 100 % espacial hacia el eje Z
beta_deg = 0;       % Giro del plano (2-3) respecto del plano (1-2) (grados)

offset_deg = 20;    % Desfase de montaje entre las horquillas del eje intermedio (grados)
                    % Homocinética si delta1 = delta2 y offset = beta

phi1_vista_deg = 30; % Posición angular del motor para dibujar las horquillas (grados)

% REFERENCIA: phi1 = 0 cuando la horquilla de entrada (eje 1) está
% contenida en el plano formado por los ejes 1 y 2 (plano XY).

%% CONVERSIÓN A RADIANES
delta1      = deg2rad(delta1_deg);
delta2      = deg2rad(delta2_deg);
beta_azimut = deg2rad(beta_deg);
offset      = deg2rad(offset_deg);
phi1_vista  = deg2rad(phi1_vista_deg);

phi1_deg = linspace(0, 360, 1000);
phi1     = deg2rad(phi1_deg);

%% CÁLCULO CINEMÁTICO
% Junta 1
num2    = cos(delta1);
den2    = 1 - sin(delta1).^2 .* cos(phi1).^2;
omega_2 = omega_1 .* (num2 ./ den2);

phi2 = atan2(sin(phi1), cos(phi1) .* cos(delta1));   % tan(phi2) = tan(phi1)/cos(delta1)
phi2 = unwrap(phi2);

% Junta 2: ángulo de la horquilla motriz (eje 2) medido desde el plano (2-3)
% (+90° porque la horquilla del eje 2 es perpendicular a la del eje 1 en la cruceta)
phi2_efectivo = phi2 + pi/2 + offset - beta_azimut;

num3    = cos(delta2);
den3    = 1 - sin(delta2).^2 .* cos(phi2_efectivo).^2;
omega_3 = omega_2 .* (num3 ./ den3);

% Rizado de la velocidad de salida
rizado = (max(omega_3) - min(omega_3)) / omega_1 * 100;
fprintf('omega_3: min = %.2f rad/s | max = %.2f rad/s | rizado = %.2f %%\n', ...
        min(omega_3), max(omega_3), rizado);

%% GEOMETRÍA 3D (BASE DE VECTORES LOCALES)
L1 = 3; L2 = 4; L3 = 3;

P0 = [0; 0; 0];
P1 = [L1; 0; 0];

% Base vectorial local del Eje 2
ux = [cos(delta1); sin(delta1); 0];    % Dirección del eje 2
uy = [-sin(delta1); cos(delta1); 0];
uz = [0; 0; 1];

P2 = P1 + L2 * ux;

% Dirección del Eje 3 (generada desde el cono y el azimut)
a3 = cos(delta2)*ux + sin(delta2)*(cos(beta_azimut)*uy + sin(beta_azimut)*uz);
P3 = P2 + L3 * a3;

%% ORIENTACIÓN DE LAS HORQUILLAS (para phi1 = phi1_vista)
a1 = [1; 0; 0];
a2 = ux;

% Rotación de Rodrigues de un vector v alrededor del eje unitario k un ángulo th
rotar = @(v, k, th) v*cos(th) + cross(k, v)*sin(th) + k*dot(k, v)*(1 - cos(th));

f1  = rotar([0; 1; 0], a1, phi1_vista);            % Horquilla de entrada (eje 1)
f2a = cross(a2, f1);  f2a = f2a / norm(f2a);       % Horquilla del eje 2 en la junta 1
f2b = rotar(f2a, a2, offset);                      % Horquilla del eje 2 en la junta 2 (con desfase)
f3  = cross(a3, f2b); f3 = f3 / norm(f3);          % Horquilla de salida (eje 3)

rad = 0.6;   % Semilongitud de los brazos de las horquillas (para el dibujo)

%% GRÁFICA 1: VISTA 3D
figure('Name', 'Esquema 3D Espacial', 'Position', [50, 100, 900, 550]);

plot3([P0(1) P1(1)], [P0(2) P1(2)], [P0(3) P1(3)], 'k-', 'LineWidth', 6, 'DisplayName', 'Eje 1 (Motor)'); hold on;
plot3([P1(1) P2(1)], [P1(2) P2(2)], [P1(3) P2(3)], 'b-', 'LineWidth', 5, 'DisplayName', 'Eje 2 (Intermedio)');
plot3([P2(1) P3(1)], [P2(2) P3(2)], [P2(3) P3(3)], 'r-', 'LineWidth', 6, 'DisplayName', 'Eje 3 (Salida)');
plot3(P1(1), P1(2), P1(3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
plot3(P2(1), P2(2), P2(3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');

% Brazos de las horquillas (segmentos que muestran su orientación real)
dibujar_horquilla = @(P, f, color, nombre) plot3([P(1)-rad*f(1) P(1)+rad*f(1)], ...
                                                 [P(2)-rad*f(2) P(2)+rad*f(2)], ...
                                                 [P(3)-rad*f(3) P(3)+rad*f(3)], ...
                                                 '-', 'Color', color, 'LineWidth', 3, 'DisplayName', nombre);
dibujar_horquilla(P1, f1,  [0.2 0.2 0.2], 'Horquilla eje 1');
dibujar_horquilla(P1, f2a, [0 0 0.8],     'Horquilla eje 2 (junta 1)');
dibujar_horquilla(P2, f2b, [0.3 0.5 1],   'Horquilla eje 2 (junta 2)');
dibujar_horquilla(P2, f3,  [0.8 0 0],     'Horquilla eje 3');

grid on; axis equal;
xlabel('Eje X'); ylabel('Eje Y'); zlabel('Eje Z');
title(sprintf('Vista 3D: \\delta_1=%.1f°, \\delta_2=%.1f°, Azimut=%.1f°, Offset=%.1f°, \\phi_1=%.1f°', ...
      delta1_deg, delta2_deg, beta_deg, offset_deg, phi1_vista_deg));
legend('Location', 'best');
view(45, 20);

%% GRÁFICA 2: VISTA PLANTA XY
figure('Name', 'Proyección en el Plano XY', 'Position', [980, 100, 750, 400]);
plot([P0(1) P1(1)], [P0(2) P1(2)], 'k-', 'LineWidth', 6); hold on;
plot([P1(1) P2(1)], [P1(2) P2(2)], 'b-', 'LineWidth', 5);
plot([P2(1) P3(1)], [P2(2) P3(2)], 'r-', 'LineWidth', 6);
plot(P1(1), P1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
plot(P2(1), P2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
grid on; axis equal; xlabel('Eje X'); ylabel('Eje Y');
title('Plano XY (Planta)');

%% GRÁFICA 3: VISTA ALZADO XZ
figure('Name', 'Proyección en el Plano XZ', 'Position', [980, 580, 750, 400]);
plot([P0(1) P1(1)], [P0(3) P1(3)], 'k-', 'LineWidth', 6); hold on;
plot([P1(1) P2(1)], [P1(3) P2(3)], 'b-', 'LineWidth', 5);
plot([P2(1) P3(1)], [P2(3) P3(3)], 'r-', 'LineWidth', 6);
plot(P1(1), P1(3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
plot(P2(1), P2(3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
grid on; axis equal; xlabel('Eje X'); ylabel('Eje Z');
title('Plano XZ (Elevación)');

%% GRÁFICA 4: CINEMÁTICA
figure('Name', 'Comportamiento de Velocidades Angulares', 'Position', [50, 700, 900, 400]);
plot(phi1_deg, omega_1 * ones(size(phi1_deg)), 'k--', 'LineWidth', 1.5); hold on;
plot(phi1_deg, omega_2, 'b-', 'LineWidth', 2);
plot(phi1_deg, omega_3, 'r-', 'LineWidth', 2);
xline(phi1_vista_deg, ':', 'Posición dibujada', 'HandleVisibility', 'off');
grid on; xlabel('Ángulo motor \phi_1 (grados)'); ylabel('Velocidad (rad/s)');
title(sprintf('Cinemática (Azimut=%.1f°, Offset=%.1f°, Rizado \\omega_3=%.2f %%)', ...
      beta_deg, offset_deg, rizado));
legend('\omega_1 (Motor)', '\omega_2 (Intermedio)', '\omega_3 (Salida)', 'Location', 'best');
xlim([0 360]);
