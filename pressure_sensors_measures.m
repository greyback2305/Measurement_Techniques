%% Static pressure vs height analysis
clear; clc; close all;

files = dir('data_press_static_*cm.csv');
N = numel(files);

h_cm  = zeros(N,1);
p_mean = zeros(N,1);
p_std  = zeros(N,1);
p_all  = cell(N,1);

for k = 1:N
    % parse height from filename: data_press_static_<H>cm.csv
    tok = regexp(files(k).name, 'static_(\d+)cm', 'tokens', 'once');
    h_cm(k) = str2double(tok{1});

    T = readtable(files(k).name);
    p = T.pressure_Pa;
    p_all{k}  = p;
    p_mean(k) = mean(p);
    p_std(k)  = std(p);
end

% sort by height
[h_cm, idx] = sort(h_cm);
p_mean = p_mean(idx);
p_std  = p_std(idx);
p_all  = p_all(idx);

h_m = h_cm/100;

%% Linear fit: p = a*h + b   -> slope a = rho*g (sign depends on convention)
[pcoef, S] = polyfit(h_m, p_mean, 1);
[p_fit, delta] = polyval(pcoef, h_m, S);
slope = pcoef(1);
intercept = pcoef(2);

g = 9.81;
rho_est = abs(slope)/g;

% R^2
ss_res = sum((p_mean - p_fit).^2);
ss_tot = sum((p_mean - mean(p_mean)).^2);
R2 = 1 - ss_res/ss_tot;

fprintf('Slope     : %.2f Pa/m\n', slope);
fprintf('Intercept : %.2f Pa\n', intercept);
fprintf('R^2       : %.4f\n', R2);
fprintf('rho est.  : %.1f kg/m^3 (assuming g = 9.81)\n', rho_est);

%% Plots
figure('Color','w','Position',[100 100 1100 750]);

% (1) raw time series, one subplot per height
subplot(2,2,1); hold on; grid on;
cmap = turbo(N);
for k = 1:N
    plot(p_all{k}, 'Color', cmap(k,:), 'DisplayName', sprintf('%d cm', h_cm(k)));
end
xlabel('sample #'); ylabel('pressure [Pa]');
title('Raw time series'); legend('Location','best');

% (2) boxplot per height
subplot(2,2,2);
G = []; P = [];
for k = 1:N
    P = [P; p_all{k}];
    G = [G; repmat(h_cm(k), numel(p_all{k}), 1)];
end
boxplot(P, G); grid on;
xlabel('height [cm]'); ylabel('pressure [Pa]');
title('Distribution per height');

% (3) mean +/- std with linear fit
subplot(2,2,[3 4]); hold on; grid on;
errorbar(h_cm, p_mean, p_std, 'o', 'MarkerFaceColor',[0.2 0.4 0.8], ...
    'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName','data (mean \pm std)');
hh = linspace(min(h_m), max(h_m), 100);
plot(hh*100, polyval(pcoef, hh), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('fit: %.1f Pa/m, R^2=%.3f', slope, R2));
xlabel('height [cm]'); ylabel('pressure [Pa]');
title(sprintf('p(h) — \\rho \\approx %.0f kg/m^3', rho_est));
legend('Location','best');