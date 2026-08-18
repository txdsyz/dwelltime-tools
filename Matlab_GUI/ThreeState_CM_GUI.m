function ThreeState_CM_GUI()
% ThreeState_CM_GUI — Three-State Markov Process Error Analysis (Simulation-based)
% How to run: Type "ThreeState_CM_GUI" directly into the MATLAB command window and press Enter.

%   Because pro = 0.5 forces pi proportional to τ, the three states share the
%   same event rate, so a single count value applies to A, B and C alike:
%       T = count * (tau_a + tau_b + tau_c)      N_pts = T * f0

    %% Main Window
    fig = uifigure('Name', 'Three-State CM Markov Simulation Analysis', ...
                   'Position', [200 100 500 760], 'Resize', 'off', ...
                   'Color', [0.97 0.97 0.97]);

    uilabel(fig, 'Text', 'CM Three-State Markov Process Uncertainty Analysis', ...
            'Position', [20 700 460 35], 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'FontColor', [0.15 0.25 0.55]);

    uipanel(fig, 'Position', [20 690 460 2], 'BackgroundColor', [0.7 0.7 0.8], 'BorderType', 'none');

    %% Input Area
    lblW = 230; fldX = 270; fldW = 140; rowH = 28; startY = 620; gap = 35;

    % Time Constants
    uilabel(fig, 'Text', '\tau_a (s):', 'Interpreter', 'tex', 'Position', [30 startY lblW rowH], 'FontSize', 12);
    fld_taua = uieditfield(fig, 'numeric', 'Position', [fldX startY+2 fldW 24], 'Value', 10);

    uilabel(fig, 'Text', '\tau_b (s):', 'Interpreter', 'tex', 'Position', [30 startY-gap lblW rowH], 'FontSize', 12);
    fld_taub = uieditfield(fig, 'numeric', 'Position', [fldX startY-gap+2 fldW 24], 'Value', 5);

    uilabel(fig, 'Text', '\tau_c (s):', 'Interpreter', 'tex', 'Position', [30 startY-2*gap lblW rowH], 'FontSize', 12);
    fld_tauc = uieditfield(fig, 'numeric', 'Position', [fldX startY-2*gap+2 fldW 24], 'Value', 2);

    uilabel(fig, 'Text', 'Sampling frequency f_0 (Hz):', 'Interpreter', 'tex', ...
            'Position', [30 startY-3*gap lblW rowH], 'FontSize', 11);
    fld_f0 = uieditfield(fig, 'numeric', 'Position', [fldX startY-3*gap+2 fldW 24], 'Value', 1);

    uilabel(fig, 'Text', 'Monte Carlo Iterations (M):', 'Interpreter', 'tex', ...
            'Position', [30 startY-4*gap lblW rowH], 'FontSize', 11);
    fld_M = uieditfield(fig, 'numeric', 'Position', [fldX startY-4*gap+2 fldW 24], 'Value', 50);

    %% Input mode selector
    uilabel(fig, 'Text', 'Experiment scale — choose ONE input:', ...
        'Position', [30 startY-5.7*gap 440 rowH], 'FontSize', 11, ...
        'FontWeight', 'bold', 'FontColor', [0.15 0.25 0.55]);

    bg = uibuttongroup(fig, ...
        'Position',        [30 startY-6.5*gap 440 28], ...
        'BorderType',      'none', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'SelectionChangedFcn', @(~,~) updateInputMode());

    rb_T     = uiradiobutton(bg, 'Text', 'Time T (s)',        'Position', [0   3 100 22], 'FontSize', 11);
    rb_count = uiradiobutton(bg, 'Text', 'Dwell count',       'Position', [290 3 150 22], 'FontSize', 11);

    tipCount = ['Expected (true) number of complete dwells per state. Because ' ...
                'p = 0.5 the three states share the same event rate, so one value ' ...
                'applies to A, B and C. Note the counts reported below are the ' ...
                'events actually recovered by Viterbi, which are fewer.'];

    % Time T field
    lbl_T = uilabel(fig, 'Text', 'Time T (s):', 'Interpreter', 'tex', ...
        'Position', [30 startY-7.5*gap lblW rowH], 'FontSize', 11);
    fld_T = uieditfield(fig, 'numeric', 'Position', [fldX startY-7.5*gap+2 fldW 24], ...
        'Value', 500, 'Limits', [eps Inf]);


    % Dwell count field
    lbl_count = uilabel(fig, 'Text', 'Expected dwell average count:', ...
        'Position', [30 startY-7.5*gap lblW rowH], 'FontSize', 11, 'Tooltip', tipCount);
    fld_count = uieditfield(fig, 'numeric', 'Position', [fldX startY-7.5*gap+2 fldW 24], ...
        'Value', 29.4, 'Limits', [eps Inf], 'Tooltip', tipCount);

    % default: Time T mode (matches the original GUI)
    rb_T.Value = true;
    updateInputMode();

    function updateInputMode()
        lbl_T.Visible     = onoff(rb_T.Value);      fld_T.Visible     = onoff(rb_T.Value);
        lbl_count.Visible = onoff(rb_count.Value);  fld_count.Visible = onoff(rb_count.Value);
    end

    %% Calculate Button
    btn = uibutton(fig, 'Text', '▶  Calculate', 'Position', [165 300 170 40], ...
                   'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.15 0.35 0.75], ...
                   'FontColor', 'white', 'ButtonPushedFcn', @(~,~) runCalc());

    %% Result Panel
    resPanel = uipanel(fig, ...
        'Title',           'Results', ...
        'Position',        [20 15 460 275], ...
        'FontSize',        11, ...
        'FontWeight',      'bold', ...
        'BackgroundColor', [0.94 0.97 1.0], ...
        'BorderType',      'line');

    lbl_result = uilabel(resPanel, 'Text', '— Click "Calculate" to start —', ...
                         'Position', [15 10 430 235], 'FontSize', 11, 'WordWrap', 'on', ...
                         'Interpreter', 'tex', ...
                         'VerticalAlignment', 'top');

    %% Core Computing Logic
    function runCalc()
        btn.Enable = 'off';
        lbl_result.Text = sprintf('Simulating %d HMM trajectories... Please wait.', fld_M.Value);
        lbl_result.FontColor = [0.4 0.4 0.4];
        drawnow;

        try

            tau_a0 = fld_taua.Value;
            tau_b0 = fld_taub.Value;
            tau_c0 = fld_tauc.Value;
            f0     = fld_f0.Value;
            M      = fld_M.Value;

            % Resolve the experiment scale to a single T
            tau_sum = tau_a0 + tau_b0 + tau_c0;
            if rb_T.Value
                T = fld_T.Value;
            else
                T = fld_count.Value * tau_sum;
            end

            count_exp = T / tau_sum;      % expected complete dwells per state
            N_pts     = T * f0;

            [Rebias, Rerr, BiasA, BiasB, BiasC, StdA, StdB, StdC, Na, Nb, Nc] = ...
                run3StateSim(tau_a0, tau_b0, tau_c0, T, f0, M);


            resStr = { ...
                sprintf('Total Relative Bias: %.4f %%', Rebias), ...
                sprintf('Total Relative Variation: %.4f %%', Rerr), ...
                '-----------------------------------------------------------------', ...
                sprintf('\\tau_a Bias : %.6e ', BiasA), ...
                sprintf('\\tau_a Variation : %.6e', StdA), ...
                sprintf('\\tau_b Bias : %.6e ', BiasB), ...
                sprintf('\\tau_b Variation : %.6e ', StdB), ...
                sprintf('\\tau_c Bias : %.6e ', BiasC), ...
                sprintf('\\tau_c Variation : %.6e ', StdC), ...
             '-----------------------------------------------------------------', ...
                sprintf('Dwell A average count  : %.1f', Na), ...
                sprintf('Dwell B average count  : %.1f', Nb), ...
                sprintf('Dwell C average count  : %.1f', Nc), ...   

            };

            lbl_result.Text = resStr;
            lbl_result.FontColor = [0.05 0.35 0.05];

        catch ME
            lbl_result.Text = ['Error: ' ME.message];
            lbl_result.FontColor = [0.7 0.05 0.05];
        end
        btn.Enable = 'on';
    end
end

function s = onoff(tf)
    if tf, s = 'on'; else, s = 'off'; end
end

%%3-State Simulation Core
function [Rebi_total, Restd_total, biasa, biasb, biasc, stda, stdb, stdc, Na, Nb, Nc] = run3StateSim(tau_a0, tau_b0, tau_c0, T, f0, M)
    alpha = 0.05;
    z_alpha = -norminv(alpha/2);
    dt = 1e-1;
    N = fix(T/dt);
    sigma0 = 0;
    dtt = 1/f0;
    Ra = fix(dtt/dt);
    NN = fix(T/dtt);
    indx = (0:NN-1)*Ra+1;

    sta_num = 10;
    sta_level = (1:sta_num)*((2-0)/sta_num);
    pro = 0.5;      % symmetric exit split,fixed model assumption


    a0 = 1/tau_a0; a1 = a0*pro;
    b0 = 1/tau_b0; b1 = b0*pro;
    c0 = 1/tau_c0; c1 = c0*pro;

    Q = [-a0, a1, a0-a1; b1, -b0, b0-b1; c1, c0-c1, -c0];
    CQ = Q .* dt;
    TRANS = expm(CQ);
    EMIS = eye(3);

    A_GUESS = [1/3 1/3 1/3; 1/3 1/3 1/3; 1/3 1/3 1/3];
    B_GUESS = [1 0 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1; 0 0 1]';

    Na_sum = 0;
    Nb_sum = 0;
    Nc_sum = 0;

    for k = 1:M
        [~, states] = hmmgenerate(N, TRANS, EMIS);
        states = [states, states(M)];
        states(1) = [];
        RawI = states + normrnd(0, sigma0, [1,N]) - 1;

        Sa_RawI = RawI(indx);
        seq = zeros(1, NN);


        for i = 1:NN
            diffs = abs(Sa_RawI(i) - sta_level);
            [~, seq(i)] = min(diffs);
        end


        [A_EST, B_EST] = hmmtrain(seq, A_GUESS, B_GUESS);
        STATES = hmmviterbi(seq, A_EST, B_EST) - 1;

        % Measure Events
        swpt = 1;
        stmark = [];
        for i = 2:NN
            if STATES(i) ~= STATES(i-1)
                swpt = [swpt, i];
                stmark = [stmark, STATES(i)];
            end
        end

        Events = diff(swpt) * dtt;
        No = length(Events);

        if No <= 1
            Tau_a(k) = -1; Tau_b(k) = -1; Tau_c(k) = -1;
        elseif No == 2
            Tau_a(k) = -1; Tau_b(k) = -1; Tau_c(k) = -1;
            if stmark(1) == 0, Tau_a(k) = Events(2); end
            if stmark(1) == 1, Tau_b(k) = Events(2); end
            if stmark(1) == 2, Tau_c(k) = Events(2); end
        else
            Events(1) = [];
            stmark(end) = [];

            Event_a = Events(stmark == 0);
            Event_b = Events(stmark == 1);
            Event_c = Events(stmark == 2);

            Na_sum = Na_sum + length(Event_a);
            Nb_sum = Nb_sum + length(Event_b);
            Nc_sum = Nc_sum + length(Event_c);

            Tau_a(k) = ifempty(mean(Event_a), -1);
            Tau_b(k) = ifempty(mean(Event_b), -1);
            Tau_c(k) = ifempty(mean(Event_c), -1);

        end

    end

    Na = Na_sum / M;
    Nb = Nb_sum / M;
    Nc = Nc_sum / M;

    Tau_a(Tau_a <= 0) = [];
    Tau_b(Tau_b <= 0) = [];
    Tau_c(Tau_c <= 0) = [];


    [biasa, stda, Rebia, Restda] = calcStats(Tau_a, tau_a0, z_alpha);

    [biasb, stdb, Rebib, Restdb] = calcStats(Tau_b, tau_b0, z_alpha);

    [biasc, stdc, Rebic, Restdc] = calcStats(Tau_c, tau_c0, z_alpha);

    Restd_total = nansum([Restda, Restdb, Restdc]);
    Rebi_total  = nansum([abs(Rebia), abs(Rebib), abs(Rebic)]);
end

function val = ifempty(val, default)
    if isempty(val) || isnan(val), val = default; end
end

function [bias, std_val, Rebi, Restd] = calcStats(Tau_data, tau_true, z_alpha)
    if isempty(Tau_data)
        bias = NaN; std_val = NaN; Rebi = NaN; Restd = NaN;
    else
        tau_mean = mean(Tau_data);
        std_val = std(Tau_data);
        bias = tau_mean - tau_true;
        Restd = z_alpha * std_val / tau_true * 100;
        Rebi = bias / tau_true * 100;
    end
end
