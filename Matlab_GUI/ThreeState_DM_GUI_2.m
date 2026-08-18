function ThreeState_DM_GUI()
% ThreeState_DM_GUI — 3-State Markov Process Error Analysis (Dye-Cycling Mode)
% How to run: Type "ThreeState_DM_GUI" directly into the MATLAB command window and press Enter.

% The experiment scale can be specified in two ways: Time T (s), or dwell
% average count.  The bright duty cycle is T_on/(T_on+T_off), so
% T = count * (tau_a+tau_b+tau_c) * (T_on+T_off)/T_on
%  pro = 0.5, which makesm the three states share the same event rate, so one count value applies to all.

    %% Main Window
    fig = uifigure('Name', 'Three-State DM Mode Error Analysis', ...
                   'Position', [200 50 500 780], 'Resize', 'off', ...
                   'Color', [0.97 0.97 0.97]);

    uilabel(fig, 'Text', 'DM Three-State Markov Process Uncertainty Analysis', ...
            'Position', [20 740 460 35], 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'FontColor', [0.15 0.25 0.55]);

    uipanel(fig, 'Position', [20 730 460 2], 'BackgroundColor', [0.7 0.7 0.8], 'BorderType', 'none');

    %% Input Area
    lblW = 230; fldX = 270; fldW = 140; rowH = 28; startY = 680; gap = 35;

    % Time Constants (Measured)
    uilabel(fig, 'Text', 'Measured \tau_a (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY lblW rowH], 'FontSize', 12);
    fld_taua = uieditfield(fig, 'numeric', 'Position', [fldX startY+2 fldW 24], 'Value', 0.1);

    uilabel(fig, 'Text', 'Measured \tau_b (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-gap lblW rowH], 'FontSize', 12);
    fld_taub = uieditfield(fig, 'numeric', 'Position', [fldX startY-gap+2 fldW 24], 'Value', 0.1);

    uilabel(fig, 'Text', 'Measured \tau_c (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-2*gap lblW rowH], 'FontSize', 12);
    fld_tauc = uieditfield(fig, 'numeric', 'Position', [fldX startY-2*gap+2 fldW 24], 'Value', 1);

    % DM Specific Parameters
    uilabel(fig, 'Text', 'Time T_{on} (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-3*gap lblW rowH], 'FontSize', 11);
    fld_Ton = uieditfield(fig, 'numeric', 'Position', [fldX startY-3*gap+2 fldW 24], 'Value', 42);

    uilabel(fig, 'Text', 'Time T_{off} (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-4*gap lblW rowH], 'FontSize', 11);
    fld_Toff = uieditfield(fig, 'numeric', 'Position', [fldX startY-4*gap+2 fldW 24], 'Value', 4);

    % Experimental Parameters
    uilabel(fig, 'Text', 'Sampling freq f_0 (Hz):', 'Interpreter', 'tex', ...
            'Position', [30 startY-5*gap lblW rowH], 'FontSize', 11);
    fld_f0 = uieditfield(fig, 'numeric', 'Position', [fldX startY-5*gap+2 fldW 24], 'Value', 1);

    uilabel(fig, 'Text', 'Monte Carlo Iterations (M):', 'Interpreter', 'tex', ...
            'Position', [30 startY-6*gap lblW rowH], 'FontSize', 11);
    fld_M = uieditfield(fig, 'numeric', 'Position', [fldX startY-6*gap+2 fldW 24], 'Value', 50);

    %% Input mode selector
    uilabel(fig, 'Text', 'Experiment scale — choose ONE input:', ...
        'Position', [30 startY-7.2*gap 440 rowH], 'FontSize', 11, ...
        'FontWeight', 'bold', 'FontColor', [0.15 0.25 0.55]);

    bg = uibuttongroup(fig, ...
        'Position',        [30 startY-8*gap 440 28], ...
        'BorderType',      'none', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'SelectionChangedFcn', @(~,~) updateInputMode());

    rb_T     = uiradiobutton(bg, 'Text', 'Time T (s)',  'Position', [0   3 200 22], 'FontSize', 11);
    rb_count = uiradiobutton(bg, 'Text', 'Dwell count', 'Position', [220 3 220 22], 'FontSize', 11);

    tipCount = 'Expected number of complete dwells per state (same for A, B and C).';

    % Time T field
    lbl_T = uilabel(fig, 'Text', 'Time T (s):', 'Interpreter', 'tex', ...
        'Position', [30 startY-9*gap lblW rowH], 'FontSize', 11);
    fld_T = uieditfield(fig, 'numeric', 'Position', [fldX startY-9*gap+2 fldW 24], ...
        'Value', 11000, 'Limits', [eps Inf]);

    % Dwell count field
    lbl_count = uilabel(fig, 'Text', 'Dwell average count:', ...
        'Position', [30 startY-9*gap lblW rowH], 'FontSize', 11, 'Tooltip', tipCount);
    fld_count = uieditfield(fig, 'numeric', 'Position', [fldX startY-9*gap+2 fldW 24], ...
        'Value', 8369.6, 'Limits', [eps Inf], 'Tooltip', tipCount);

    % default: Time T mode (matches the original GUI)
    rb_T.Value = true;
    updateInputMode();

    function updateInputMode()
        lbl_T.Visible     = onoff(rb_T.Value);      fld_T.Visible     = onoff(rb_T.Value);
        lbl_count.Visible = onoff(rb_count.Value);  fld_count.Visible = onoff(rb_count.Value);
    end

    %% Calculate Button
    btn = uibutton(fig, 'Text', '▶  Calculate', 'Position', [165 290 170 40], ...
                   'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.15 0.35 0.75], ...
                   'FontColor', 'white', 'ButtonPushedFcn', @(~,~) runCalc());

    %% Result Panel
    resPanel = uipanel(fig, ...
        'Title',           'Result', ...
        'Position',        [20 15 460 265], ...
        'FontSize',        11, ...
        'FontWeight',      'bold', ...
        'BackgroundColor', [0.94 0.97 1.0], ...
        'BorderType',      'line');

    lbl_result = uilabel(resPanel, 'Text', '— Click "Calculate" to start —', ...
                         'Position', [15 10 430 225], 'FontSize', 11, 'WordWrap', 'on', ...
                         'Interpreter', 'tex', ...
                         'VerticalAlignment', 'top');

    %% Core Computing Logic
    function runCalc()
        btn.Enable = 'off';
        lbl_result.Text = 'Running...';
        lbl_result.FontColor = [0.4 0.4 0.4];
        drawnow;

        try

            tau_m_a = fld_taua.Value;
            tau_m_b = fld_taub.Value;
            tau_m_c = fld_tauc.Value;
            T_on    = fld_Ton.Value;
            T_off   = fld_Toff.Value;
            f0      = fld_f0.Value;
            M       = fld_M.Value;

            % Resolve the experiment scale to a single T
            tau_sum = tau_m_a + tau_m_b + tau_m_c;
            duty    = T_on / (T_on + T_off);          % bright fraction
            if rb_T.Value
                T = fld_T.Value;
            else
                T = fld_count.Value * tau_sum / duty;
            end

            N_pts = T * f0;

            [Rebias, Rerr, BiasA, BiasB, BiasC, StdA, StdB, StdC, Na, Nb, Nc] = ...
                eval3StateDMError(tau_m_a, tau_m_b, tau_m_c, T_on, T_off, T, f0, M);

            Rerr_str = formatOutput(Rerr, '%%');

            resStr = { ...
                sprintf('Total Expected Relative Bias: %.4f %%', Rebias), ...
                ['Total Expected Variation: ', Rerr_str], ...
                '-----------------------------------------------------------------', ...
                sprintf('\\tau_a Bias : %.6e ', BiasA), ...
                sprintf('\\tau_a Variation : %.6e', StdA), ...
                sprintf('\\tau_b Bias : %.6e ', BiasB), ...
                sprintf('\\tau_b Variation : %.6e ', StdB), ...
                sprintf('\\tau_c Bias : %.6e ', BiasC), ...
                sprintf('\\tau_c Variation : %.6e ', StdC), ...
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

function str = formatOutput(val, unit)
    if isinf(val) || isnan(val)
        str = 'N/A (Insufficient events)';
    else
        if strcmp(unit, '%%')
            str = sprintf('%.4f %%', val);
        else
            str = sprintf('%.6e %s', val, unit);
        end
    end
end

%%3-State DM 
function [Rebi_total, Restd_total, biasa, biasb, biasc, stda, stdb, stdc, Na, Nb, Nc] = eval3StateDMError(tau_a, tau_b, tau_c, T_on, T_off, T, f0, M)
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
    sta_level = (1:sta_num)*((3-0)/sta_num);
    pro = 0.5;   
    
    
    TB_win = (1/T_on)/(1/T_on+1/T_off);
    TA_win = (1/T_off)/(1/T_on+1/T_off);
    TE_win = exp(-(1/T_on+1/T_off)*dt);
    T_TRANS = [TB_win+TA_win*TE_win, TA_win-TA_win*TE_win; ...
               TB_win-TB_win*TE_win, TA_win+TB_win*TE_win];
    T_EMIS = [1, 0; 0, 1];
    

    a0 = 1/tau_a; a1 = a0*pro;
    b0 = 1/tau_b; b1 = b0*pro;
    c0 = 1/tau_c; c1 = c0*pro;
    Q = [-a0, a1, a0-a1; b1, -b0, b0-b1; c1, c0-c1, -c0];
    CQ = Q .* dt;
    TRANS = expm(CQ);
    EMIS = eye(3);
    
    Tau_a_sim = zeros(1, M);
    Tau_b_sim = zeros(1, M);
    Tau_c_sim = zeros(1, M);
    

    A_GUESS = ones(4,4)/4;
    B_GUESS = [1 0 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0; 0 1 0 0; ...
               0 0 1 0; 0 0 1 0; 0 0 1 0; 0 0 0 1; 0 0 0 1]';
    Na_sum = 0;
    Nb_sum = 0;
    Nc_sum = 0;

    for k = 1:M
     
        [~, states] = hmmgenerate(N, TRANS, EMIS);
        states = [states, states(end)];
        states(1) = [];
        
        
        [~, Tstates] = hmmgenerate(N, T_TRANS, T_EMIS);
        Tstates = [Tstates, Tstates(end)];
        Tstates(1) = [];
        WindowI = Tstates - 1;
        
    
        RawI = states .* WindowI + normrnd(0,sigma0,[1,N]);
        

        Sa_RawI = RawI(indx);
    
        [~, seq] = min(abs(Sa_RawI(:) - sta_level), [], 2);
        seq = seq'; 
   
        [A_EST, B_EST] = hmmtrain(seq, A_GUESS, B_GUESS);
        STATES = hmmviterbi(seq, A_EST, B_EST) - 1;
        
    
        swpt = 1;
        Inista = STATES(1);
        for i = 2:NN
            if STATES(i) ~= STATES(i-1)
               swpt = [swpt, i];
               Inista = [Inista, STATES(i)];
            end
        end
        
        if isempty(Inista)
            Tau_a_sim(k) = NaN; Tau_b_sim(k) = NaN; Tau_c_sim(k) = NaN;
            continue;
        end
        
        if Inista(1) == 0
            blank = [1, 2];
        elseif length(Inista) >= 2 && Inista(2) == 0
            blank = [1, 2, 3];
        else
            blank = 1;
        end
        
        for i = 3:length(swpt)
            if Inista(i) == 0
               blank = [blank, i-1, i, i+1];
            end
        end
        
        Events = diff(swpt) * dtt;
        Inista(end) = [];
        
        blank = unique(blank(blank <= length(Events)));
        
        Events(blank) = [];
        Inista(blank) = [];
        
   
        Event_a = Events(Inista == 1);
        Event_b = Events(Inista == 2);
        Event_c = Events(Inista == 3);

        Na_sum = Na_sum + length(Event_a); 
        Nb_sum = Nb_sum + length(Event_b);
        Nc_sum = Nc_sum + length(Event_c);
            
        Tau_a_sim(k) = ifempty(mean(Event_a), NaN);
        Tau_b_sim(k) = ifempty(mean(Event_b), NaN);
        Tau_c_sim(k) = ifempty(mean(Event_c), NaN);
    end
    Na = Na_sum / M;
    Nb = Nb_sum / M;
    Nc = Nc_sum / M;

    Tau_a_sim(isnan(Tau_a_sim) | Tau_a_sim <= 0) = [];
    Tau_b_sim(isnan(Tau_b_sim) | Tau_b_sim <= 0) = [];
    Tau_c_sim(isnan(Tau_c_sim) | Tau_c_sim <= 0) = [];
    
    [biasa, stda, Rebia, Restda] = calcStats(Tau_a_sim, tau_a, z_alpha);
    [biasb, stdb, Rebib, Restdb] = calcStats(Tau_b_sim, tau_b, z_alpha);
    [biasc, stdc, Rebic, Restdc] = calcStats(Tau_c_sim, tau_c, z_alpha);
    
    Restd_total = nansum([Restda, Restdb, Restdc]);
    Rebi_total  = nansum([abs(Rebia), abs(Rebib), abs(Rebic)]);
end

function val = ifempty(val, default)
    if isempty(val) || isnan(val), val = default; end
end

function [bias, std_val, Rebi, Restd] = calcStats(Tau_data, tau_meas, z_alpha)
    if isempty(Tau_data)
        bias = NaN; std_val = NaN; Rebi = NaN; Restd = NaN;
    else
        tau_sim_mean = mean(Tau_data);
        std_val = std(Tau_data);
        bias = tau_sim_mean - tau_meas;
        Restd = z_alpha * std_val / tau_meas * 100;
        Rebi = bias / tau_meas * 100;
    end
end