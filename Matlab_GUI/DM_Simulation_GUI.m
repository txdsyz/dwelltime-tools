function DM_simulation_GUI()
% DM_Sim_GUI — Two-State Markov Process, Monte Carlo simulation (DyeCycling Mode)
% How to run: Type "DM_simulation_GUI" directly into the MATLAB command window and press Enter.


    %% Main Window
    fig = uifigure('Name', 'DM Two-State Markov —  Simulation', ...
                   'Position', [200 60 500 730], 'Resize', 'off', ...
                   'Color', [0.97 0.97 0.97]);

    uilabel(fig, 'Text', 'DM Two-State Markov Process Simulation', ...
            'Position', [20 680 460 26], 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'FontColor', [0.15 0.25 0.55]);

    uipanel(fig, 'Position', [20 648 460 2], 'BackgroundColor', [0.7 0.7 0.8], 'BorderType', 'none');

    %% Input Area
    lblW = 230; fldX = 270; fldW = 140; rowH = 28; startY = 600; gap = 35;

    uilabel(fig, 'Text', '\tau_a (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY lblW rowH], 'FontSize', 12);
    fld_taua = uieditfield(fig, 'numeric', 'Position', [fldX startY+2 fldW 24], ...
            'Value', 10, 'Limits', [eps Inf]);

    uilabel(fig, 'Text', '\tau_b (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-gap lblW rowH], 'FontSize', 12);
    fld_taub = uieditfield(fig, 'numeric', 'Position', [fldX startY-gap+2 fldW 24], ...
            'Value', 10, 'Limits', [eps Inf]);

    tipOn = 'Mean duration of a on window .';
    uilabel(fig, 'Text', 'T_{on} (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-2*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipOn);
    fld_Ton = uieditfield(fig, 'numeric', 'Position', [fldX startY-2*gap+2 fldW 24], ...
            'Value', 50, 'Limits', [eps Inf], 'Tooltip', tipOn);

    tipOff = 'Mean duration of a off window.';
    uilabel(fig, 'Text', 'T_{off} (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-3*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipOff);
    fld_Toff = uieditfield(fig, 'numeric', 'Position', [fldX startY-3*gap+2 fldW 24], ...
            'Value', 5, 'Limits', [eps Inf], 'Tooltip', tipOff);

    tipFull = ['Total length of one trace. Number of bright windows per trace is ' ...
               'about T_full/(T_on+T_off). '];
    uilabel(fig, 'Text', 'Trace length T_{full} (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-4*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipFull);
    fld_Tfull = uieditfield(fig, 'numeric', 'Position', [fldX startY-4*gap+2 fldW 24], ...
            'Value', 1100, 'Limits', [eps Inf], 'Tooltip', tipFull, ...
            'ValueChangedFcn', @(~,~) updateEstimate());

    uilabel(fig, 'Text', 'Sampling frequency f_0 (Hz):', 'Interpreter', 'tex', ...
            'Position', [30 startY-5*gap lblW rowH], 'FontSize', 12);
    fld_f0 = uieditfield(fig, 'numeric', 'Position', [fldX startY-5*gap+2 fldW 24], ...
            'Value', 1, 'Limits', [eps Inf], 'ValueChangedFcn', @(~,~) updateEstimate());

    tipN = ['Number of traces per simulation run.'];
    uilabel(fig, 'Text', 'Number of traces N:', 'Interpreter', 'tex', ...
            'Position', [30 startY-6*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipN);
    fld_N = uieditfield(fig, 'numeric', 'Position', [fldX startY-6*gap+2 fldW 24], ...
            'Value', 2, 'Limits', [1 Inf], 'RoundFractionalValues', 'on', ...
            'Tooltip', tipN, 'ValueChangedFcn', @(~,~) updateEstimate());

    tipM = ['Number of independent simulation runs.'];
    uilabel(fig, 'Text', 'Monte Carlo runs M:', 'Interpreter', 'tex', ...
            'Position', [30 startY-7*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipM);
    fld_M = uieditfield(fig, 'numeric', 'Position', [fldX startY-7*gap+2 fldW 24], ...
            'Value', 5, 'Limits', [1 Inf], 'RoundFractionalValues', 'on', ...
            'Tooltip', tipM, 'ValueChangedFcn', @(~,~) updateEstimate());

    %% Cost estimate
    lbl_est = uilabel(fig, 'Text', '', ...
        'Position', [30 328 440 20], 'FontSize', 10, ...
        'FontColor', [0.45 0.45 0.45], 'FontAngle', 'italic');
    updateEstimate();

    function updateEstimate()
        dt    = 1e-3;
        N_tr  = round(fld_Tfull.Value/dt);
        total = 2 * N_tr * fld_N.Value * fld_M.Value;
        nwin  = fld_Tfull.Value / (fld_Ton.Value + fld_Toff.Value);

        if total > 2e8
            lbl_est.FontColor = [0.7 0.05 0.05];
        else
            lbl_est.FontColor = [0.45 0.45 0.45];
        end
    end

    %% Calculate Button
    btn = uibutton(fig, 'Text', '▶  Run simulation', 'Position', [155 270 190 38], ...
                   'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.15 0.35 0.75], ...
                   'FontColor', 'white', 'ButtonPushedFcn', @(~,~) runCalc());

    %% Result Panel
    resPanel = uipanel(fig, ...
        'Title',           'Simulation Results', ...
        'Position',        [20 15 460 245], ...
        'FontSize',        11, ...
        'FontWeight',      'bold', ...
        'BackgroundColor', [0.94 0.97 1.0], ...
        'BorderType',      'line');

    lbl_result = uilabel(resPanel, 'Text', '— Click "Run simulation" to start —', ...
                         'Position', [15 10 430 205], 'FontSize', 11, 'WordWrap', 'on', ...
                         'Interpreter', 'tex', ...
                         'VerticalAlignment', 'top');

    %% Core Computing Logic
    function runCalc()
        btn.Enable = 'off';
        lbl_result.Text = 'Simulating... Please wait.';
        lbl_result.FontColor = [0.4 0.4 0.4];
        drawnow;

        try
            ta0    = fld_taua.Value;   tb0   = fld_taub.Value;
            T_on   = fld_Ton.Value;    T_off = fld_Toff.Value;
            T_full = fld_Tfull.Value;  f0    = fld_f0.Value;
            N      = fld_N.Value;      M     = fld_M.Value;

            % constants
            alpha    = 0.05;  z_alpha = -norminv(alpha/2);
            dt       = 1e-3;
            dtt      = 1/f0;
            Ra       = round(dtt/dt);
            minValid = 5;
            sta_num  = 10;  sta_level = (1:sta_num)*(2/sta_num);   % 0.2 .. 2.0
            A_GUESS  = ones(3)/3;
            B_GUESS  = [1 0 0;1 0 0;1 0 0; 0 1 0;0 1 0;0 1 0; 0 0 1;0 0 1;0 0 1;0 0 1]';
            EMIS2    = [1 0; 0 1];
            T_EMIS   = [1 0; 0 1];

            if Ra < 1
                error('Sampling interval 1/f0 is shorter than the simulation step dt = 1e-3 s. Reduce f0 below 1000 Hz.');
            end

            % on/off window transition matrix (fixed, independent of tau)
            TBw = (1/T_on)/(1/T_on+1/T_off);  TAw = (1/T_off)/(1/T_on+1/T_off);
            TEw = exp(-(1/T_on+1/T_off)*dt);
            T_TRANS = [TBw+TAw*TEw TAw-TAw*TEw; TBw-TBw*TEw TAw+TBw*TEw];

            [ta, tb, sa, sb, ba, bb, Rsa, Rsb, Rba, Rbb, nva, nEvt] = ...
                runDMSim(ta0, tb0, N, M, T_full, dt, Ra, dtt, z_alpha, ...
                         sta_level, A_GUESS, B_GUESS, EMIS2, T_TRANS, T_EMIS, ...
                         minValid);

            if isnan(ta)
                resStr = { ...
                    sprintf('Not enough valid runs: %d of %d (minimum %d).', nva, M, minValid), ...
                    'After discarding dwells adjacent to dark segments, too few events', ...
                    'survive. Try a longer T_{full}, a larger T_{on}, or a higher f_0.' };
            else
                resStr = { ...
                    sprintf('Total Relative Bias: %.4f %%', abs(Rba) + abs(Rbb)), ...
                    sprintf('Total Relative Variation: %.4f %%', Rsa + Rsb), ...
                    '-----------------------------------------------------------------', ...
                    sprintf('\\tau_a  true : %.4g s   measured : %.4g s', ta0, ta), ...
                    sprintf('\\tau_a  Bias : %.4e    Variation : %.4e', ba, sa), ...
                    sprintf('\\tau_b  true : %.4g s   measured : %.4g s', tb0, tb), ...
                    sprintf('\\tau_b  Bias : %.4e    Variation : %.4e', bb, sb), ...
 };
            end

            lbl_result.Text = resStr;
            lbl_result.FontColor = [0.05 0.35 0.05];

        catch ME
            lbl_result.Text = ['Error: ' ME.message];
            lbl_result.FontColor = [0.7 0.05 0.05];
        end
        btn.Enable = 'on';
    end
end

%% Single-point simulation 
function [ta,tb,sa,sb,ba,bb,Rsa,Rsb,Rba,Rbb,nva,nEvt] = ...
         runDMSim(ta0, tb0, N, M, T_full, dt, Ra, dtt, z_alpha, ...
                  sta_level, A_GUESS, B_GUESS, EMIS2, T_TRANS, T_EMIS, ...
                  minValid)

    ta=nan;tb=nan;sa=nan;sb=nan;ba=nan;bb=nan;Rsa=nan;Rsb=nan;Rba=nan;Rbb=nan;
    nva=0; nEvt=0;

    alph = 1/ta0;  beta = 1/tb0;
    Bc = beta/(alph+beta); Ac = alph/(alph+beta); Ec = exp(-(alph+beta)*dt);
    TRANS = [Bc+Ac*Ec Ac-Ac*Ec; Bc-Bc*Ec Ac+Bc*Ec];

    N_tr = round(T_full/dt);

    Tau_a = nan(1,M);  Tau_b = nan(1,M);  cntA = nan(1,M);

    for k = 1:M
        allEa = [];  allEb = [];
        for tr = 1:N
            [~,states]  = hmmgenerate(N_tr, TRANS, EMIS2);
            states = [states(2:end), states(end)];
            [~,Tstates] = hmmgenerate(N_tr, T_TRANS, T_EMIS);
            WindowI = Tstates - 1;            % 0=off, 1=on
            RawI = states .* WindowI;         % 0=off, 1=a, 2=b
            states=[]; Tstates=[]; WindowI=[];

            idx = 1:Ra:N_tr;  Sa = RawI(idx);  RawI=[];
            NN = numel(Sa);
            seq = zeros(1,NN);
            for i = 1:NN, [~,seq(i)] = min(abs(Sa(i)-sta_level)); end
 
            try
                [A_EST,B_EST] = hmmtrain(seq, A_GUESS, B_GUESS);
                STATES = hmmviterbi(seq, A_EST, B_EST) - 1;   % 0=off,1=a,2=b
            catch
                continue;
            end
         
            swpt = 1;  Inista = STATES(1);
            for i = 2:NN
                if STATES(i) ~= STATES(i-1)
                    swpt = [swpt, i];  Inista = [Inista, STATES(i)];
                end
            end
            if numel(swpt) < 3, continue; end
            Events = diff(swpt)*dtt;      % length = numel(swpt)-1
            Inista(end) = [];             %  Events
            % remove off dwell and before & after dwell
            % remove the first dwell
            blank = 1;
            if Inista(1) == 0
                blank = [blank, 2];
            elseif numel(Inista) >= 2 && Inista(2) == 0
                blank = [blank, 2, 3];
            end
            for i = 3:numel(Inista)
                if Inista(i) == 0, blank = [blank, i-1, i, i+1]; end
            end
            blank = unique(blank);
            blank = blank(blank >= 1 & blank <= numel(Events));
            Events(blank) = [];  Inista(blank) = [];
            if isempty(Events), continue; end
            %  a / b 
            Ea = Events(Inista == 1);
            Eb = Events(Inista == 2);
            allEa = [allEa, Ea];  allEb = [allEb, Eb];
        end
        if ~isempty(allEa), Tau_a(k) = mean(allEa); end
        if ~isempty(allEb), Tau_b(k) = mean(allEb); end
        cntA(k) = numel(allEa);
    end

    nva  = sum(~isnan(Tau_a));
    nEvt = mean(cntA,'omitnan');
    if nva >= minValid
        ta=mean(Tau_a,'omitnan'); tb=mean(Tau_b,'omitnan');
        sa=std(Tau_a,'omitnan');  sb=std(Tau_b,'omitnan');
        ba=ta-ta0;  bb=tb-tb0;
        Rsa=z_alpha*sa/ta0*100;  Rsb=z_alpha*sb/tb0*100;
        Rba=ba/ta0*100;          Rbb=bb/tb0*100;
    end
end