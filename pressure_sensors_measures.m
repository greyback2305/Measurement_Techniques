%% Static pressure vs height analysis
clear; clc; close all;

%% Load and parse files
files = dir('data_press_static_*cm.csv');
N = numel(files);
h_cm   = zeros(N,1);
p_mean = zeros(N,1);
p_std  = zeros(N,1);
p_all  = cell(N,1);

for k = 1:N
    tok = regexp(files(k).name, 'static_(\d+p\d+|\d+)cm', 'tokens', 'once');
    h_cm(k) = str2double(strrep(tok{1}, 'p', '.'));
    T = readtable(files(k).name);
    p = T.pressure_Pa;
    p_all{k}  = p;
    p_mean(k) = mean(p);
    p_std(k)  = std(p);
end

%% Sort by height
[h_cm, idx] = sort(h_cm);
p_mean = p_mean(idx);
p_std  = p_std(idx);
p_all  = p_all(idx);
h_mm = h_cm * 10;
h_m  = h_cm / 100;

%% Linear fit (h in metres)
[pcoef, S]  = polyfit(h_m, p_mean, 1);
[p_fit, ~]  = polyval(pcoef, h_m, S);
slope       = pcoef(1);
intercept   = pcoef(2);
g           = 9.81;
rho_est     = abs(slope) / g;
rho_water   = 1000;

ss_res = sum((p_mean - p_fit).^2);
ss_tot = sum((p_mean - mean(p_mean)).^2);
R2 = 1 - ss_res / ss_tot;

fprintf('Slope     : %.2f Pa/m\n',  slope);
fprintf('Intercept : %.2f Pa\n',    intercept);
fprintf('R^2       : %.4f\n',       R2);
fprintf('rho est.  : %.1f kg/m^3\n', rho_est);
fprintf('rho error : %.1f%%\n', abs(rho_est - rho_water)/rho_water * 100);

%% Box half-width in mm (30% of min spacing)
bw = 0.30 * min(diff(h_mm));

cmap = turbo(N);

%% ---- Figure -----------------------------------------------
figure('Color','w','Position',[80 80 1200 500]);

%% -- Subplot 1 : Raw time series -----------------------------------------
subplot(1,2,1); hold on; grid on; box on;

for k = 1:N
    plot(p_all{k}, 'Color', cmap(k,:), ...
        'DisplayName', sprintf('%.1f mm', h_mm(k)));
end
xlabel('Sample #');
ylabel('Pressure [Pa]');
title('Raw pressure time series');
legend('Location','best','NumColumns',2);

%% -- Subplot 2 : Manual boxplot at real h_mm + fit -----------------------
subplot(1,2,2); hold on; grid on; box on;

for k = 1:N
    p  = p_all{k};
    q  = quantile(p, [0 0.25 0.50 0.75 1]);   % min, Q1, median, Q3, max
    x  = h_mm(k);

    % Whiskers (min–max)
    plot([x x], [q(1) q(5)], 'k-', 'LineWidth', 1);
    % Whisker caps
    plot([x-bw/2 x+bw/2], [q(1) q(1)], 'k-', 'LineWidth', 1);
    plot([x-bw/2 x+bw/2], [q(5) q(5)], 'k-', 'LineWidth', 1);
    % IQR box
    patch([x-bw x+bw x+bw x-bw], [q(2) q(2) q(4) q(4)], ...
        'w', 'EdgeColor','k', 'LineWidth', 1.2);
    % Median line
    plot([x-bw x+bw], [q(3) q(3)], 'k-', 'LineWidth', 2);
end


% Fit and theory lines in real mm / m space
hh_mm = linspace(min(h_mm), max(h_mm), 200);
hh_m  = hh_mm / 1000;

plot(hh_mm, polyval(pcoef, hh_m), 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('fit: %.0f Pa/m  (\\rho=%.0f kg/m^3, R^2=%.4f)', ...
    slope, rho_est, R2));

p_theory = rho_water * g * hh_m + intercept;
plot(hh_mm, p_theory, 'b--', 'LineWidth', 2, ...
    'DisplayName', sprintf('theory: %.0f Pa/m  (\\rho_w=%.0f kg/m^3)', ...
    rho_water*g, rho_water));

% Mean markers
plot(h_mm, p_mean, 'o', ...
    'MarkerFaceColor', [0.15 0.35 0.75], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, 'LineWidth', 1.2, ...
    'DisplayName', 'mean');

% Dummy patch for box legend entry

% Fit and theory lines in real mm / m space
hh_mm = linspace(min(h_mm), max(h_mm), 200);
hh_m  = hh_mm / 1000;

h1 = plot(hh_mm, polyval(pcoef, hh_m), 'r-', 'LineWidth', 2);

p_theory = rho_water * g * hh_m + intercept;
h2 = plot(hh_mm, p_theory, 'b--', 'LineWidth', 2);

h3 = plot(h_mm, p_mean, 'o', ...
    'MarkerFaceColor', [0.15 0.35 0.75], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, 'LineWidth', 1.2);

legend([h1 h2 h3], ...
    {sprintf('fit: %.0f Pa/m  (\\rho=%.0f kg/m^3, R^2=%.4f)', slope, rho_est, R2), ...
     sprintf('theory: %.0f Pa/m  (\\rho_w=%.0f kg/m^3)', rho_water*g, rho_water), ...
     'mean'}, ...
    'Location','best');


xlabel('Height [mm]');
ylabel('Pressure [Pa]');
title(sprintf('p(h) — \\rho_{est} = %.0f vs \\rho_w = %.0f kg/m^3', rho_est, rho_water));
legend('Location','best');

sgtitle('Static pressure vs height', 'FontSize', 13, 'FontWeight', 'bold');