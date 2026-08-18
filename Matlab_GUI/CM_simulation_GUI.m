function CM_simulation_GUI()
% How to run: Type "CM_simulation_GUI" directly into the MATLAB command window and press Enter.

    %% Main Window
    fig = uifigure('Name', 'CM Two-State Markov — Simulation', ...
                   'Position', [200 80 500 700], 'Resize', 'off', ...
                   'Color', [0.97 0.97 0.97]);

    uilabel(fig, 'Text', 'CM Two-State Markov Process Simulation', ...
            'Position', [20 650 460 26], 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'FontColor', [0.15 0.25 0.55]);

    uipanel(fig, 'Position', [20 618 460 2], 'BackgroundColor', [0.7 0.7 0.8], 'BorderType', 'none');

    %% Input Area
    lblW = 230; fldX = 270; fldW = 140; rowH = 28; startY = 570; gap = 35;

    uilabel(fig, 'Text', '\tau_a (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY lblW rowH], 'FontSize', 12);
    fld_taua = uieditfield(fig, 'numeric', 'Position', [fldX startY+2 fldW 24], ...
            'Value', 10, 'Limits', [eps Inf]);

    uilabel(fig, 'Text', '\tau_b (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-gap lblW rowH], 'FontSize', 12);
    fld_taub = uieditfield(fig, 'numeric', 'Position', [fldX startY-gap+2 fldW 24], ...
            'Value', 10, 'Limits', [eps Inf]);

    tipL = ['Mean length of one trace. Traces are terminated at random, ' ...
            'so individual lengths are drawn from an exponential ' ...
            'distribution with this mean. ' ...
            'intervals are discarded.'];
    uilabel(fig, 'Text', 'Mean trace length L (s):', 'Interpreter', 'tex', ...
            'Position', [30 startY-2*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipL);
    fld_L = uieditfield(fig, 'numeric', 'Position', [fldX startY-2*gap+2 fldW 24], ...
            'Value', 100, 'Limits', [eps Inf], 'Tooltip', tipL);

    uilabel(fig, 'Text', 'Sampling frequency f_0 (Hz):', 'Interpreter', 'tex', ...
            'Position', [30 startY-3*gap lblW rowH], 'FontSize', 12);
    fld_f0 = uieditfield(fig, 'numeric', 'Position', [fldX startY-3*gap+2 fldW 24], ...
            'Value', 1, 'Limits', [eps Inf]);

    tipN = ['Number of traces per simulation run.'];
    uilabel(fig, 'Text', 'Number of traces N:', 'Interpreter', 'tex', ...
            'Position', [30 startY-4*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipN);
    fld_N = uieditfield(fig, 'numeric', 'Position', [fldX startY-4*gap+2 fldW 24], ...
            'Value', 5, 'Limits', [1 Inf], 'RoundFractionalValues', 'on', 'Tooltip', tipN);

    tipM = ['Number of independent simulation runs. Bias and spread are computed ' ...
            'across these runs. The production grid uses 50; a smaller value runs '];
    uilabel(fig, 'Text', 'Monte Carlo runs M:', 'Interpreter', 'tex', ...
            'Position', [30 startY-5*gap lblW rowH], 'FontSize', 12, 'Tooltip', tipM);
    fld_M = uieditfield(fig, 'numeric', 'Position', [fldX startY-5*gap+2 fldW 24], ...
            'Value', 10, 'Limits', [1 Inf], 'RoundFractionalValues', 'on', 'Tooltip', tipM);

    %% Calculate Button
    btn = uibutton(fig, 'Text', '▶  Run simulation', 'Position', [155 330 190 38], ...
                   'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.15 0.35 0.75], ...
                   'FontColor', 'white', 'ButtonPushedFcn', @(~,~) runCalc());

    %% Result Panel
    resPanel = uipanel(fig, ...
        'Title',           'Simulation Results', ...
        'Position',        [20 15 460 300], ...
        'FontSize',        11, ...
        'FontWeight',      'bold', ...
        'BackgroundColor', [0.94 0.97 1.0], ...
        'BorderType',      'line');

    lbl_result = uilabel(resPanel, 'Text', '— Click "Run simulation" to start —', ...
                         'Position', [15 10 430 260], 'FontSize', 11, 'WordWrap', 'on', ...
                         'Interpreter', 'tex', ...
                         'VerticalAlignment', 'top');

    %% Core Computing Logic
    function runCalc()
        btn.Enable = 'off';
        lbl_result.Text = 'Simulating... Please wait.';
        lbl_result.FontColor = [0.4 0.4 0.4];
        drawnow;

        try
            ta0 = fld_taua.Value;  tb0 = fld_taub.Value;
            L   = fld_L.Value;     f0  = fld_f0.Value;
            N   = fld_N.Value;     M   = fld_M.Value;

            % ---- fixed constants, matching CM_2D_prod ----
            alpha    = 0.05;  z_alpha = -norminv(alpha/2);
            dt       = 1e-3;
            dtt      = 1/f0;
            Ra       = round(dtt/dt);
            sigma0   = 0;
            minValid = 5;
            sta_num  = 10;  sta_level = (1:sta_num)*(1/sta_num);
            A_GUESS  = [1/2 1/2; 1/2 1/2];
            B_GUESS  = [1 0;1 0;1 0;1 0;1 0; 0 1;0 1;0 1;0 1;0 1]';
            EMIS     = [1 0; 0 1];

            if Ra < 1
                error('Sampling interval 1/f0 is shorter than the simulation step dt = 1e-3 s. Reduce f0 below 1000 Hz.');
            end

            [ta, tb, sa, sb, ba, bb, Rsa, Rsb, Rba, Rbb, nva, nEvt] = ...
                runCMSim(ta0, tb0, N, M, L, dt, Ra, dtt, z_alpha, sigma0, ...
                         sta_level, A_GUESS, B_GUESS, EMIS, minValid);

            if isnan(ta)
                resStr = { ...
                    sprintf('Not enough valid runs: %d of %d (minimum %d).', nva, M, minValid), ...
                    'The traces are too short, or the dwells are unresolvable at this f_0.', ...
                    'Try a longer L, a higher f_0, or larger N.' };
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

%% ============ Single-point simulation (algorithm from oneCell) ============
function [ta,tb,sa,sb,ba,bb,Rsa,Rsb,Rba,Rbb,nva,nEvt] = ...
         runCMSim(ta0, tb0, N, M, L, dt, Ra, dtt, z_alpha, sigma0, ...
                  sta_level, A_GUESS, B_GUESS, EMIS, minValid)

    ta=nan;tb=nan;sa=nan;sb=nan;ba=nan;bb=nan;Rsa=nan;Rsb=nan;Rba=nan;Rbb=nan;
    nva=0; nEvt=0;

    alph = 1/ta0;  beta = 1/tb0;
    Bc = beta/(alph+beta); Ac = alph/(alph+beta); Ec = exp(-(alph+beta)*dt);
    TRANS = [Bc+Ac*Ec Ac-Ac*Ec; Bc-Bc*Ec Ac+Bc*Ec];

    Tau_a = nan(1,M);  Tau_b = nan(1,M);  cntA = nan(1,M);

    for k = 1:M
        % 事件汇总:把 N 条 trace 的所有 dwell 混成一大池,取一次 mean
        allEa = [];  allEb = [];
        for tr = 1:N
            Lon = exprnd(L);  N_tr = round(Lon/dt);
            if N_tr < 20*Ra, continue; end
            [~,states] = hmmgenerate(N_tr, TRANS, EMIS);
            states = [states(2:end), states(end)];
            RawI = states + normrnd(0,sigma0,[1,N_tr]) - 1;
            idx = 1:Ra:N_tr;  Sa = RawI(idx);  NN = numel(Sa);
            seq = zeros(1,NN);
            for i = 1:NN, [~,seq(i)] = min(abs(Sa(i)-sta_level)); end
            [A_EST,B_EST] = hmmtrain(seq, A_GUESS, B_GUESS);
            STATES = hmmviterbi(seq, A_EST, B_EST) - 1;
            swpt = 1;
            for i = 2:NN
                if STATES(i) ~= STATES(i-1), swpt = [swpt, i]; end
            end
            if numel(swpt) < 3, continue; end
            Events = diff(swpt)*dtt;  Events(1) = [];
            No = numel(Events);
            if mod(No,2)==0, ia=(0:No/2-1)*2+1; ib=(1:No/2)*2;
            else, ia=(0:floor(No/2))*2+1; ib=(1:floor(No/2))*2; end
            if STATES(swpt(2))==0, Ea=Events(ia); Eb=Events(ib);
            else, Eb=Events(ia); Ea=Events(ib); end
            allEa = [allEa, Ea];  allEb = [allEb, Eb];
        end
        % 一个仿真实例 k:所有 trace 的事件汇总后取一次 mean
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