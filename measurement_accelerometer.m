clear; clc; close all;

% --- Load data ---
filename = 'data_acc_n1_n1_off_1.csv';

if isfile(filename)
    raw = readlines(filename);
    data = [];
    for i = 1:length(raw)
        line = strtrim(raw(i));
        if line == "", continue; end
        idx = strfind(line, ' -> ');
        if ~isempty(idx)
            line = extractAfter(line, idx(1) + 3);
        end
        vals = str2double(strsplit(line, ','));
        if length(vals) == 7 && ~any(isnan(vals))
            data = [data; vals];
        end
    end
    fprintf('Loaded %d samples from %s\n', size(data,1), filename);
else
    error('data.csv not found.');
end





t     = data(:,1);
Ax    = data(:,2);
Ay    = data(:,3);
Az    = data(:,4);
Aabs  = data(:,5);
roll  = data(:,6);
pitch = data(:,7);


% =========================================================
% --- Axis calibration from static measurement ---
% =========================================================
cal_file = 'data_acc_zero.csv';

if isfile(cal_file)
    raw_cal = readlines(cal_file);
    cal_data = [];
    for i = 1:length(raw_cal)
        line = strtrim(raw_cal(i));
        if line == "", continue; end
        idx = strfind(line, ' -> ');
        if ~isempty(idx)
            line = extractAfter(line, idx(1) + 3);
        end
        vals = str2double(strsplit(line, ','));
        if length(vals) == 7 && ~any(isnan(vals))
            cal_data = [cal_data; vals];
        end
    end

    % Remove outliers from calibration data
    cal_Ax   = cal_data(:,2);
    cal_Ay   = cal_data(:,3);
    cal_Az   = cal_data(:,4);

    bad_cal = abs(zscore(cal_Ax)) > 3 | abs(zscore(cal_Ay)) > 3 | abs(zscore(cal_Az)) > 3;
    cal_Ax(bad_cal) = []; cal_Ay(bad_cal) = []; cal_Az(bad_cal) = [];

    % Static mean = bias on each axis
    % One axis should be ~1000 mg (gravity), others ~0
    bias_Ax = mean(cal_Ax);
    bias_Ay = mean(cal_Ay);
    bias_Az = mean(cal_Az);

    % Identify gravity axis
    [~, grav_idx] = max(abs([bias_Ax, bias_Ay, bias_Az]));
    names = {'Ax', 'Ay', 'Az'};
    fprintf('\n=== Calibration ===\n');
    fprintf('Static bias — Ax: %.1f mg, Ay: %.1f mg, Az: %.1f mg\n', bias_Ax, bias_Ay, bias_Az);
    biases = [bias_Ax, bias_Ay, bias_Az];
    fprintf('Gravity axis: %s (%.1f mg)\n', names{grav_idx}, biases(grav_idx));

    % Subtract bias, keep 1000 mg on gravity axis
    gravity = 1000;
    offsets = [bias_Ax, bias_Ay, bias_Az];
    offsets(grav_idx) = offsets(grav_idx) - sign(offsets(grav_idx)) * gravity;

    Ax = Ax - offsets(1);
    Ay = Ay - offsets(2);
    Az = Az - offsets(3);

    fprintf('Applied offsets — Ax: %.1f mg, Ay: %.1f mg, Az: %.1f mg\n', offsets(1), offsets(2), offsets(3));
    fprintf('Post-calibration means — Ax: %.1f, Ay: %.1f, Az: %.1f mg\n', mean(Ax), mean(Ay), mean(Az));


        % --- Remap axes so gravity is always +Z ---
    ax_vals = [Ax, Ay, Az];
    
    % Reorder so gravity axis becomes Az
    Ax = ax_vals(:, mod(grav_idx,   3) + 1);   % first non-gravity axis
    Ay = ax_vals(:, mod(grav_idx+1, 3) + 1);   % second non-gravity axis
    Az = ax_vals(:, grav_idx);                  % gravity axis → Z
    
    % Flip Z to be positive if needed
    if mean(Az) < 0
        Az = -Az;
    end

fprintf('Remapped axes — gravity is now +Az\n');
fprintf('Post-remap means — Ax: %.1f, Ay: %.1f, Az: %.1f mg\n', mean(Ax), mean(Ay), mean(Az));
else
    warning('Calibration file not found — skipping axis correction.');
end


roll  = atan2d(Ay, Az);
pitch = atan2d(-Ax, sqrt(Ay.^2 + Az.^2));

% Shift roll wrapping point away from ±180
roll = roll - mean(roll);
roll  = unwrap(roll  * pi/180) * 180/pi;
roll = mod((roll+180), 360)-180 ;

pitch = pitch - mean(pitch);
pitch  = unwrap(pitch  * pi/180) * 180/pi;
pitch = mod((pitch+180), 360)-180 ;

% --- Remove bad samples ---
bad = (Ax == 0) | (Az == 0) | abs(zscore(Aabs)) > 3;
fprintf('Removed %d bad samples\n', sum(bad));
t(bad)=[]; Ax(bad)=[]; Ay(bad)=[]; Az(bad)=[];
Aabs(bad)=[]; roll(bad)=[]; pitch(bad)=[];

% --- Sampling frequency ---
dt = mean(diff(t));
fs = 1 / dt;
fprintf('Sampling rate: %.1f Hz\n', fs);

% --- Resample to uniform time grid (fixes jitter) ---
t_uniform = (t(1):dt:t(end))';
Ax    = interp1(t, Ax,    t_uniform, 'linear');
Ay    = interp1(t, Ay,    t_uniform, 'linear');
Az    = interp1(t, Az,    t_uniform, 'linear');
Aabs  = interp1(t, Aabs,  t_uniform, 'linear');
roll  = interp1(t, roll,  t_uniform, 'linear');
pitch = interp1(t, pitch, t_uniform, 'linear');
t = t_uniform;

% --- Statistics ---
fprintf('\n=== Statistics (mg) ===\n');
fprintf('%6s %8s %8s %8s %8s\n', '', 'Ax', 'Ay', 'Az', '|A|');
fprintf('%6s %8.2f %8.2f %8.2f %8.2f\n', 'Mean', mean(Ax), mean(Ay), mean(Az), mean(Aabs));
fprintf('%6s %8.2f %8.2f %8.2f %8.2f\n', 'Std',  std(Ax),  std(Ay),  std(Az),  std(Aabs));
fprintf('%6s %8.2f %8.2f %8.2f %8.2f\n', 'Min',  min(Ax),  min(Ay),  min(Az),  min(Aabs));
fprintf('%6s %8.2f %8.2f %8.2f %8.2f\n', 'Max',  max(Ax),  max(Ay),  max(Az),  max(Aabs));

% =========================================================
% --- FFT ---
% =========================================================
signals   = {Ax,   Ay,   Az,   roll,   pitch};
sig_names = {'Ax', 'Ay', 'Az', 'Roll', 'Pitch'};
sig_units = {'mg', 'mg', 'mg', 'deg',  'deg'};
f_min     = 0.5;   % ignore drift below this frequency

figure('Name', 'FFT Analysis', 'Position', [100 100 1000 700]);
fprintf('\n=== Dominant frequencies ===\n');
for k = 1:5
    sig = signals{k} - mean(signals{k});   % remove DC
    N   = length(sig);
    win = hanning(N);                       % Hanning window
    Y   = abs(fft(sig .* win)) / sum(win);
    Y   = Y(1:floor(N/2)+1);
    Y(2:end-1) = 2 * Y(2:end-1);
    f   = fs * (0:floor(N/2)) / N;

    % Find dominant frequency above f_min
    Y_valid = Y;
    Y_valid(f < f_min) = 0;
    [~, idx] = max(Y_valid);
    fprintf('Dominant frequency in %s: %.2f Hz\n', sig_names{k}, f(idx));

    subplot(3, 2, k);
    plot(f, Y, 'LineWidth', 1.2);
    xline(f(idx), 'r--', sprintf('%.2f Hz', f(idx)), 'LabelVerticalAlignment', 'bottom');
    xlabel('Frequency (Hz)'); ylabel(['Amplitude (' sig_units{k} ')']);
    title(['FFT — ' sig_names{k}]); grid on;
    xlim([0 fs/2]);
end

% =========================================================
% --- Filter ---
% =========================================================
cutoff_hz = 5;   % <-- tune this after checking FFT
order     = 4;
[b, a]    = butter(order, cutoff_hz / (fs/2), 'low');

Ax_f    = filtfilt(b, a, Ax);
Ay_f    = filtfilt(b, a, Ay);
Az_f    = filtfilt(b, a, Az);
roll_f  = filtfilt(b, a, roll);
pitch_f = filtfilt(b, a, pitch);

% =========================================================
% --- Raw vs Filtered ---
% =========================================================
figure('Name', 'Raw vs Filtered', 'Position', [100 100 1000 700]);

subplot(3,1,1);
plot(t, Az, 'b', 'LineWidth', 0.5); hold on;
plot(t, Az_f, 'r', 'LineWidth', 1.8);
legend('Raw Ax', 'Filtered Ax'); xlabel('Time (s)'); ylabel('mg');
title(['Ax — Low-pass ' num2str(cutoff_hz) ' Hz']); grid on;

subplot(3,1,2);
plot(t, roll, 'b', 'LineWidth', 0.5); hold on;
plot(t, roll_f, 'r', 'LineWidth', 1.8);
legend('Raw Roll', 'Filtered Roll'); xlabel('Time (s)'); ylabel('deg');
title(['Roll — Low-pass ' num2str(cutoff_hz) ' Hz']); grid on;

subplot(3,1,3);
plot(t, pitch, 'b', 'LineWidth', 0.5); hold on;
plot(t, pitch_f, 'r', 'LineWidth', 1.8);
legend('Raw Pitch', 'Filtered Pitch'); xlabel('Time (s)'); ylabel('deg');
title(['Pitch — Low-pass ' num2str(cutoff_hz) ' Hz']); grid on;

% =========================================================
% --- Position in X/Y from pitch and roll ---
% =========================================================

% Use filtered angles
% pitch → forward/back tilt (X direction)
% roll  → left/right tilt  (Y direction)
% Project onto a unit sphere: X = sin(pitch), Y = sin(roll)
X_pos = sind(pitch_f);
Y_pos = sind(roll_f);

figure('Name', 'Angular Position', 'Position', [100 100 1000 700]);

% Time series of X and Y position
subplot(3,1,1);
plot(t, X_pos, 'b', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('sin(pitch)');
title('X position over time (pitch)'); grid on;

subplot(3,1,2);
plot(t, Y_pos, 'r', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('sin(roll)');
title('Y position over time (roll)'); grid on;

subplot(3,1,3);

% Draw line colored by time using patch trick
patch([X_pos; nan], [Y_pos; nan], [t; nan], ...
    'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.2);
colorbar; xlabel('X — sin(pitch)'); ylabel('Y — sin(roll)');
title('2D angular trajectory (color = time)');
axis equal; grid on; hold on;

% Mark start and end
plot(X_pos(1),   Y_pos(1),   'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(X_pos(end), Y_pos(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
legend('', 'Start', 'End', 'Location', 'best');


% =========================================================
% --- 3D angular trajectory + time subplots ---
% =========================================================
figure('Name', 'Angular Trajectory 3D', 'Position', [100 100 1200 800]);

% 3D plot: X, Y, time
subplot(2,2,[1 3]);
patch([X_pos; nan], [Y_pos; nan], [t; nan], [t; nan], ...
    'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.2);
colorbar;
hold on;
plot3(X_pos(1),   Y_pos(1),   t(1),   'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot3(X_pos(end), Y_pos(end), t(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('X — sin(pitch)'); ylabel('Y — sin(roll)'); zlabel('Time (s)');
title('3D angular trajectory (Z = time)');
legend('', 'Start', 'End'); grid on; view(45, 30);

% X vs time
subplot(2,2,2);
plot(t, X_pos, 'b', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('sin(pitch)');
title('X (pitch) vs time'); grid on;

% Y vs time
subplot(2,2,4);
plot(t, Y_pos, 'r', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('sin(roll)');
title('Y (roll) vs time'); grid on;

% =========================================================
% --- Acceleration + Angles overview ---
% =========================================================
figure('Name', 'Accelerometer Data', 'Position', [100 100 900 700]);

subplot(3,1,1);
plot(t, Ax, '-b', t, Ay, '-g', t, Az, '-r', 'LineWidth', 1.2);
hold on; plot(t, Aabs, '--m', 'LineWidth', 1.5);
legend('Ax','Ay','Az','|A|'); xlabel('Time (s)'); ylabel('mg');
title('Acceleration components'); grid on;

subplot(3,1,2);
plot(t, roll, '-b', t, pitch, '-r', 'LineWidth', 1.2);
legend('Roll','Pitch'); xlabel('Time (s)'); ylabel('deg');
title('Tilt angles'); grid on;

subplot(3,1,3);
bar(t, zscore(Ay), 'FaceColor', [0.2 0.6 0.8]);
yline(2,'r--','+2σ'); yline(-2,'r--','-2σ');
xlabel('Time (s)'); ylabel('Z-score');
title('Ay z-score (outlier check)'); grid on;

% --- Export ---
results = table(t, Ax, Ay, Az, Aabs, roll, pitch, Ax_f, Ay_f, Az_f, roll_f, pitch_f);
writetable(results, 'results_export.csv');
fprintf('\nExported to results_export.csv\n');