%% Batch accelerometer analysis: per-file + cross-config comparison
clear; clc; close all;

%% Settings
cal_file  = 'data_acc_zero.csv';
data_glob = 'data_acc_*.csv';
out_dir   = 'figs_per_file';
cutoff_hz = 5;
order     = 4;
f_min     = 0.5;
save_per_file_png = true;     % write a PNG per file
show_per_file     = false;    % don't pop 22 figures

if ~exist(out_dir,'dir'); mkdir(out_dir); end
if ~show_per_file; set(0,'DefaultFigureVisible','off'); end

%% Calibration (once)
cal = load_acc_csv(cal_file);
[offsets, grav_idx] = compute_calibration(cal);
ax_names = {'Ax','Ay','Az'};
fprintf('Calibration: offsets=[%.1f %.1f %.1f] mg, gravity on %s\n\n', ...
    offsets(1), offsets(2), offsets(3), ax_names{grav_idx});

%% File list
all_files = dir(data_glob);
files = all_files(~strcmp({all_files.name}, cal_file));
N = numel(files);
fprintf('Found %d data files.\n', N);

%% Process every file
varnames = {'file','config','trial','n_samples','fs', ...
            'Ax_std','Ay_std','Az_std','Aabs_mean','Aabs_std', ...
            'roll_std','pitch_std', ...
            'fdom_Ax','fdom_Ay','fdom_Az','fdom_Aabs','fdom_roll','fdom_pitch'};
vt = [{'string','string','double'}, repmat({'double'},1,numel(varnames)-3)];
S = table('Size',[N numel(varnames)],'VariableTypes',vt,'VariableNames',varnames);

% keep spectra in memory for cross-config overlay
spectra_all = cell(N,1);

for k = 1:N
    fname = files(k).name;
    [cfg, trial] = parse_filename(fname);
    fprintf('[%2d/%2d] %s (config=%s, trial=%d)\n', k, N, fname, cfg, trial);

    raw = load_acc_csv(fname);
    if size(raw,1) < 50
        warning('  Too few samples, skipping.'); continue;
    end
    R = analyze_acc(raw, offsets, grav_idx, cutoff_hz, order, f_min);

    S.file(k)       = fname;
    S.config(k)     = cfg;
    S.trial(k)      = trial;
    S.n_samples(k)  = R.N;
    S.fs(k)         = R.fs;
    S.Ax_std(k)     = std(R.Ax_f);
    S.Ay_std(k)     = std(R.Ay_f);
    S.Az_std(k)     = std(R.Az_f);
    S.Aabs_mean(k)  = mean(R.Aabs);
    S.Aabs_std(k)   = std(R.Aabs);
    S.roll_std(k)   = std(R.roll_f);
    S.pitch_std(k)  = std(R.pitch_f);
    S.fdom_Ax(k)    = R.f_dom.Ax;
    S.fdom_Ay(k)    = R.f_dom.Ay;
    S.fdom_Az(k)    = R.f_dom.Az;
    S.fdom_Aabs(k)  = R.f_dom.Aabs;
    S.fdom_roll(k)  = R.f_dom.roll;
    S.fdom_pitch(k) = R.f_dom.pitch;
    spectra_all{k}  = R.spectra;

    if save_per_file_png
        plot_per_file(R, fname, cutoff_hz, ...
            fullfile(out_dir, [erase(fname,'.csv') '.png']));
    end
end

% drop empty rows from skipped files
S = S(S.n_samples > 0, :);

writetable(S, 'summary_by_file.csv');

%% Cross-config summary
set(0,'DefaultFigureVisible','on');

G = groupsummary(S, 'config', {'mean','std'}, ...
    {'Aabs_std','Ax_std','Ay_std','Az_std','roll_std','pitch_std','fdom_Aabs'});
writetable(G, 'summary_by_config.csv');
disp(G);

plot_config_bars(S);
plot_config_spectra(S, spectra_all);

fprintf('\nDone. Wrote summary_by_file.csv, summary_by_config.csv, %s/*.png\n', out_dir);

%% =========================================================
%% Local functions
%% =========================================================
function data = load_acc_csv(fname)
    if ~isfile(fname); error('File not found: %s', fname); end
    raw = readlines(fname);
    data = nan(numel(raw), 7);
    n = 0;
    for i = 1:numel(raw)
        line = strtrim(raw(i));
        if line == ""; continue; end
        idx = strfind(line, ' -> ');
        if ~isempty(idx); line = extractAfter(line, idx(1)+3); end
        v = str2double(strsplit(line, ','));
        if numel(v) == 7 && ~any(isnan(v))
            n = n+1; data(n,:) = v;
        end
    end
    data = data(1:n,:);
end

function [offsets, grav_idx] = compute_calibration(cal)
    A = cal(:,2:4);
    bad = any(abs(zscore(A)) > 3, 2);
    A(bad,:) = [];
    bias = mean(A);
    [~, grav_idx] = max(abs(bias));
    offsets = bias;
    offsets(grav_idx) = offsets(grav_idx) - sign(offsets(grav_idx)) * 1000;
end

function [cfg, trial] = parse_filename(fname)
    s = erase(erase(fname,'.csv'),'data_acc_');
    tok = regexp(s, '^(.*)_(\d+)$', 'tokens', 'once');
    if isempty(tok); cfg = string(s); trial = 0;
    else;            cfg = string(tok{1}); trial = str2double(tok{2}); end
end

function R = analyze_acc(raw, offsets, grav_idx, cutoff_hz, order, f_min)
    t = raw(:,1); A = raw(:,2:4); Aabs = raw(:,5);
    A = A - offsets;                              % bias removal
    perm = [mod(grav_idx,3)+1, mod(grav_idx+1,3)+1, grav_idx];
    A = A(:, perm);                               % gravity -> last axis
    if mean(A(:,3)) < 0; A(:,3) = -A(:,3); end
    Ax = A(:,1); Ay = A(:,2); Az = A(:,3);

    roll  = atan2d(Ay, Az);
    pitch = atan2d(-Ax, sqrt(Ay.^2 + Az.^2));

    bad = (Ax==0 & Az==0) | abs(zscore(Aabs)) > 3;
    t(bad)=[]; Ax(bad)=[]; Ay(bad)=[]; Az(bad)=[];
    Aabs(bad)=[]; roll(bad)=[]; pitch(bad)=[];

    dt = mean(diff(t)); fs = 1/dt;
    tu = (t(1):dt:t(end))';
    Ax = interp1(t,Ax,tu); Ay = interp1(t,Ay,tu); Az = interp1(t,Az,tu);
    Aabs = interp1(t,Aabs,tu);
    roll = interp1(t,roll,tu); pitch = interp1(t,pitch,tu);
    t = tu;

    [b,a] = butter(order, cutoff_hz/(fs/2), 'low');
    Ax_f    = filtfilt(b,a,Ax);
    Ay_f    = filtfilt(b,a,Ay);
    Az_f    = filtfilt(b,a,Az);
    roll_f  = filtfilt(b,a,roll);
    pitch_f = filtfilt(b,a,pitch);

    sigs = struct('Ax',Ax,'Ay',Ay,'Az',Az,'Aabs',Aabs,'roll',roll,'pitch',pitch);
    fn = fieldnames(sigs);
    f_dom = struct(); spectra = struct();
    for i = 1:numel(fn)
        s = sigs.(fn{i}) - mean(sigs.(fn{i}));
        Nlen = numel(s); win = hanning(Nlen);
        Y = abs(fft(s.*win))/sum(win);
        Y = Y(1:floor(Nlen/2)+1);
        Y(2:end-1) = 2*Y(2:end-1);
        f = fs*(0:floor(Nlen/2))/Nlen;
        Yv = Y; Yv(f<f_min) = 0;
        [~, idx] = max(Yv);
        f_dom.(fn{i}) = f(idx);
        spectra.(fn{i}) = struct('f',f(:),'Y',Y(:));
    end

    R = struct('t',t,'Ax',Ax,'Ay',Ay,'Az',Az,'Aabs',Aabs, ...
        'roll',roll,'pitch',pitch, ...
        'Ax_f',Ax_f,'Ay_f',Ay_f,'Az_f',Az_f, ...
        'roll_f',roll_f,'pitch_f',pitch_f, ...
        'fs',fs,'N',numel(t),'f_dom',f_dom,'spectra',spectra);
end

function plot_per_file(R, fname, cutoff_hz, outpath)
    fig = figure('Position',[100 100 1200 800]);
    subplot(3,2,1);
    plot(R.t,R.Ax,'b',R.t,R.Ay,'g',R.t,R.Az,'r'); hold on;
    plot(R.t,R.Aabs,'--m');
    legend('Ax','Ay','Az','|A|'); xlabel('t (s)'); ylabel('mg');
    title('Acceleration'); grid on;

    subplot(3,2,2);
    plot(R.t,R.roll,'b',R.t,R.pitch,'r');
    legend('roll','pitch'); xlabel('t (s)'); ylabel('deg');
    title('Tilt'); grid on;

    subplot(3,2,3);
    plot(R.t,R.Az,'b',R.t,R.Az_f,'r','LineWidth',1.2);
    legend('Az raw',sprintf('Az LP %g Hz',cutoff_hz));
    xlabel('t (s)'); ylabel('mg'); title('Az filtered'); grid on;

    sp = R.spectra.Aabs;
    subplot(3,2,4);
    plot(sp.f, sp.Y); xline(R.f_dom.Aabs,'r--',sprintf('%.2f Hz',R.f_dom.Aabs));
    xlim([0 R.fs/2]); xlabel('f (Hz)'); ylabel('|A|');
    title('FFT |A|'); grid on;

    sp = R.spectra.roll;
    subplot(3,2,5);
    plot(sp.f, sp.Y); xline(R.f_dom.roll,'r--',sprintf('%.2f Hz',R.f_dom.roll));
    xlim([0 R.fs/2]); xlabel('f (Hz)'); ylabel('roll');
    title('FFT roll'); grid on;

    sp = R.spectra.pitch;
    subplot(3,2,6);
    plot(sp.f, sp.Y); xline(R.f_dom.pitch,'r--',sprintf('%.2f Hz',R.f_dom.pitch));
    xlim([0 R.fs/2]); xlabel('f (Hz)'); ylabel('pitch');
    title('FFT pitch'); grid on;

    sgtitle(strrep(fname,'_','\_'));
    exportgraphics(fig, outpath, 'Resolution', 120);
    close(fig);
end

function plot_config_bars(S)
    cfgs = unique(S.config);
    nc = numel(cfgs);
    metrics = {'Aabs_std','Ax_std','Ay_std','Az_std','roll_std','pitch_std','fdom_Aabs'};
    titles  = {'std |A| [mg]','std Ax [mg]','std Ay [mg]','std Az [mg]', ...
               'std roll [deg]','std pitch [deg]','dominant f |A| [Hz]'};

    figure('Position',[100 100 1300 800],'Name','Config comparison');
    for m = 1:numel(metrics)
        subplot(3,3,m); hold on; grid on;
        mu = zeros(nc,1); sd = zeros(nc,1);
        for i = 1:nc
            v = S.(metrics{m})(S.config == cfgs(i));
            mu(i) = mean(v); sd(i) = std(v);
        end
        bar(1:nc, mu, 'FaceColor',[0.3 0.5 0.8]);
        errorbar(1:nc, mu, sd, 'k.','LineWidth',1);
        set(gca,'XTick',1:nc,'XTickLabel',cellstr(cfgs),'XTickLabelRotation',35);
        title(titles{m});
    end
    sgtitle('Cross-config metrics (mean \pm std across trials)');
end


function plot_config_spectra(S, spectra_all)
    keep  = ["n0_n0", "n1_n1", "n1_n1_off", "n1_n0"];
    cfgs  = intersect(unique(S.config), keep, 'stable');
    % reorder to tell a clean story
    [~, ord] = ismember(["n0_n0","n1_n1","n1_n1_off","n1_n0"], cfgs);
    cfgs = cfgs(ord(ord>0));
    nc    = numel(cfgs);
    cmap  = lines(nc);

    ncols = 2;
    nrows = ceil(nc / ncols);

    figure('Name', 'Vibration spectra by config', ...
           'Position', [100 100 1200 800]);

    for i = 1:nc
        idx   = find(S.config == cfgs(i));
        if isempty(idx); continue; end

        col   = cmap(i,:);
        f_ref = spectra_all{idx(1)}.Aabs.f;
        nf    = numel(f_ref);
        nt    = numel(idx);
        Y_mat = zeros(nt, nf);

        for j = 1:nt
            sp = spectra_all{idx(j)}.Aabs;
            Y_mat(j,:) = interp1(sp.f, sp.Y, f_ref, 'linear', 0);
        end

        Y_mean = mean(Y_mat, 1);
        Y_std  = std(Y_mat, 0, 1);

        subplot(nrows, ncols, i);
        hold on; grid on;

        for j = 1:nt
            plot(f_ref, Y_mat(j,:), 'Color', [col 0.35], ...
                 'LineWidth', 0.8, 'DisplayName', sprintf('trial %d', j));
        end

        fill([f_ref; flipud(f_ref)], ...
             [Y_mean + Y_std, fliplr(Y_mean - Y_std)]', ...
             col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
             'DisplayName', '\pm 1 std');

        plot(f_ref, Y_mean, 'Color', col, 'LineWidth', 2.5, ...
             'DisplayName', 'mean');

        xline(mean(S.fdom_Aabs(idx)), '--', 'Color', col*0.6, 'LineWidth', 1.6, ...
              'DisplayName', sprintf('mean f_{dom} = %.2f Hz', mean(S.fdom_Aabs(idx))));

        xlabel('Frequency (Hz)');
        ylabel('|A| [mg]');
        title(sprintf('Config: %s  (%d trials)', strrep(cfgs(i),'_','\_'), nt));
        legend('Location','northeast');
        xlim([0, max(S.fs)/2]);
    end

    sgtitle('Nozzle geometry effect on vibration (mean \pm std)');
end