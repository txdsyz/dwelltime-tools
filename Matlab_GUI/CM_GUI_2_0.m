function CM_GUI()

% CM MarkovAnalysisGUI  —  Two-State Markov Process Error Analysis
% How to run: Type "CM_GUI" directly into the MATLAB command window and press Enter.
%
% Dependency files (must be placed in the same folder as this file):
%   bM.m, bT.m, VarM.m, VarT.m
%
% NEW: The experiment scale can be specified in TWO ways (choose one):
%   (1) Total number of sampling points N_pts, or
%   (2) Dwell average count  (a single value, applied to BOTH state A and B).
% The two are mathematically interchangeable; only one is entered at a time.


    %% main window 
    fig = uifigure( ...
        'Name',     'Two-State CM Markov Process — Uncertainty Analysis', ...
        'Position', [200 120 480 700], ... 
        'Resize',   'off', ...
        'Color',    [0.97 0.97 0.97]);

    %% Title 
    uilabel(fig, ...
        'Text',                'CM Two-State Markov Process Uncertainty Analysis', ...
        'Position',            [20 650 440 35], ...
        'FontSize',            15, ...
        'FontWeight',          'bold', ...
        'HorizontalAlignment', 'center', ...
        'FontColor',           [0.15 0.25 0.55]);

    %% Divider
    uipanel(fig, 'Position', [20 640 440 2], 'BackgroundColor', [0.7 0.7 0.8], 'BorderType', 'none');

    %% Input field
    lblW = 200; fldX = 240; fldW = 140; rowH = 32; startY = 590;
    gap  = 45;
    
    % tau_a
    uilabel(fig, 'Text', '\tau_a (s):', 'Interpreter', 'tex',...
        'Position', [30 startY lblW+10 rowH], 'FontSize', 12);
    fld_taua = uieditfield(fig, 'numeric', ...
        'Position', [fldX startY+4 fldW 26], ...
        'Value', 10, 'Limits', [1e-9 Inf], 'FontSize', 11);
        
    % tau_b
    uilabel(fig, 'Text', '\tau_b (s):', 'Interpreter', 'tex', ...
        'Position', [30 startY-gap lblW+10 rowH], 'FontSize', 12);
    fld_taub = uieditfield(fig, 'numeric', ...
        'Position', [fldX startY-gap+4 fldW 26], ...
        'Value', 10, 'Limits', [1e-9 Inf], 'FontSize', 11);
        
    % f0 
    uilabel(fig, 'Text', 'Sampling frequency f₀ (Hz):', ...
        'Position', [30 startY-2*gap lblW+10 rowH], 'FontSize', 11);
    fld_f0 = uieditfield(fig, 'numeric', ...
        'Position', [fldX startY-2*gap+4 fldW 26], ...
        'Value', 100, 'Limits', [1e-9 Inf], 'FontSize', 11);

    %% ── Input mode selector ────────────────────────────────────────────────
    uilabel(fig, 'Text', 'Experiment scale — choose ONE input:', ...
        'Position', [30 420 420 rowH], 'FontSize', 11, ...
        'FontWeight', 'bold', 'FontColor', [0.15 0.25 0.55]);

    bg = uibuttongroup(fig, ...
        'Position',        [30 388 420 28], ...
        'BorderType',      'none', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'SelectionChangedFcn', @(~,~) updateInputMode());

    rb_Npts = uiradiobutton(bg, 'Text', 'Total sampling points (N)', ...
        'Position', [0 3 200 22], 'FontSize', 11);
    rb_count = uiradiobutton(bg, 'Text', 'Dwell average count', ...
        'Position', [220 3 200 22], 'FontSize', 11);

    % N_pts field
    lbl_Npts = uilabel(fig, 'Text', 'Total number of sampling points N pts:', ...
        'Position', [30 345 lblW+10 rowH], 'FontSize', 11);
    fld_Npts = uieditfield(fig, 'numeric', ...
        'Position', [fldX 349 fldW 26], ...
        'Value', 10000, 'Limits', [1 Inf], 'FontSize', 11);

    % Dwell count field (single value, applies to BOTH A and B)
    lbl_count = uilabel(fig, 'Text', 'Dwell average count (same for A & B):', ...
        'Position', [30 345 lblW+10 rowH], 'FontSize', 11);
    fld_count = uieditfield(fig, 'numeric', ...
        'Position', [fldX 349 fldW 26], ...
        'Value', 4.5, 'Limits', [1e-9 Inf], 'FontSize', 11);

    % note under the count field
    lbl_note = uilabel(fig, ...
        'Text', 'Note: one count value is used for both state A and state B.', ...
        'Position', [30 322 420 20], 'FontSize', 10, ...
        'FontColor', [0.45 0.45 0.45], 'FontAngle', 'italic');

    % default: N_pts mode
    rb_Npts.Value = true;
    updateInputMode();

    function updateInputMode()
        useNpts = rb_Npts.Value;
        if useNpts
            lbl_Npts.Visible  = 'on';  fld_Npts.Visible  = 'on';
            lbl_count.Visible = 'off'; fld_count.Visible = 'off';
            lbl_note.Visible  = 'off';
        else
            lbl_Npts.Visible  = 'off'; fld_Npts.Visible  = 'off';
            lbl_count.Visible = 'on';  fld_count.Visible = 'on';
            lbl_note.Visible  = 'on';
        end
    end

    %% Calculate button
    btn = uibutton(fig, ...
        'Text',             '▶  Calculate', ...
        'Position',         [155 265 170 38], ...
        'FontSize',         13, ...
        'FontWeight',       'bold', ...
        'BackgroundColor',  [0.15 0.35 0.75], ...
        'FontColor',        'white', ...
        'ButtonPushedFcn',  @(~,~) runCalc());

    %% Result panel 
    resPanel = uipanel(fig, ...
        'Title',           'Result', ...
        'Position',        [20 15 440 210], ...
        'FontSize',        11, ...
        'FontWeight',      'bold', ...
        'BackgroundColor', [0.94 0.97 1.0], ...
        'BorderType',      'line');
        
    lbl_result = uilabel(resPanel, ...
        'Text',                '— Click "Calculate" to view the results —', ...
        'Position',            [10 5 415 180], ...
        'FontSize',            11, ...
        'WordWrap',            'on', ...
        'Interpreter',         'tex', ...
        'VerticalAlignment',   'top', ...
        'HorizontalAlignment', 'left', ...
        'FontColor',           [0.4 0.4 0.4]);

    %% Core Computing
    function runCalc()
        % Read input from the interface
        tau_a  = fld_taua.Value;
        tau_b  = fld_taub.Value;
        f0     = fld_f0.Value;

        btn.Enable = 'off';
        lbl_result.Text      = 'Calculating......';
        lbl_result.FontColor = [0.4 0.4 0.4];
        drawnow;

        try
            % Determine T from whichever input mode is active
            if rb_Npts.Value
                N_pts = fld_Npts.Value;
                T = N_pts / f0;
            else
                count_target = fld_count.Value;
                T = T_from_count(count_target, tau_a, tau_b);
            end

            [Rebias_total, Rerr_total, Bias_a, Bias_b, Var_a, Var_b, N] = ...
                computeMarkov(tau_a, tau_b, f0, T);

            lbl_result.Interpreter = 'tex';

            % Show the derived counterpart for cross-checking
            if rb_Npts.Value
                derived_str = sprintf('Derived dwell average count : %.4f', N);
            else
                derived_str = sprintf('Derived total sampling points N : %.0f', T * f0);
            end

            resStr = { ...
                sprintf('Total Rebias: %.4f %%', Rebias_total), ...
                sprintf('Total Variation: %.4f %%', Rerr_total), ...
                '--------------------------------------------------', ...
                sprintf('\\tau_a Bias : %.6e', Bias_a), ...
                sprintf('\\tau_a Variation : %.6e', Var_a), ...
                sprintf('\\tau_b Bias : %.6e', Bias_b), ...
                sprintf('\\tau_b Variation : %.6e', Var_b), ...
                '--------------------------------------------------', ...
                derived_str, ...
            };

            lbl_result.Text = resStr;
            lbl_result.FontColor = [0.05 0.35 0.05];
        catch ME
            lbl_result.Text = ['Error: ' ME.message];
            lbl_result.FontColor = [0.7 0.05 0.05];
        end
        btn.Enable = 'on';
    end

end % ── end of GUI function ──────────────────────────────────────────────


%% ── count -> T inversion (discard_T depends on T, so iterate) ─────────────
function T = T_from_count(count_target, tau_a, tau_b)
    T = count_target * (tau_a + tau_b);   % initial guess ignoring discard
    for iter = 1:200
        discard_T = computeDiscard(T, tau_a, tau_b);
        T_new = count_target * (tau_a + tau_b) + discard_T;
        if abs(T_new - T) < 1e-13
            T = T_new;
            return;
        end
        T = T_new;
    end
end

function discard_T = computeDiscard(T, tau_a, tau_b)
    A = tau_a/(tau_a+tau_b);
    B = tau_b/(tau_a+tau_b);
    E = exp(-T*(tau_a+tau_b)/(tau_a*tau_b));
    Pa = A*(A+B*E);  Pb = A*(B-B*E);
    Pc = B*(A-A*E);  Pd = B*(B+A*E);
    discard_T = Pa*tau_a + Pd*tau_b + (Pb+Pc)*(tau_a+tau_b)/2;
end


%% ── core computation (now parametrised by T instead of N_pts) ─────────────
function [Rebias_total, Rerr_total, Bias_a, Bias_b, Std_a, Std_b, N] = ...
         computeMarkov(tau_a, tau_b, f0, T)
     
    alpha = 0.05; 
    z_alpha = -norminv(alpha/2);
    
    delta_t = 1/f0;
    cover   = 0.9;

    discard_T = computeDiscard(T, tau_a, tau_b);
    N = max(0, (T - discard_T)/(tau_a+tau_b));
    
    %%  tau_a 
    MissNb = missN(tau_b, delta_t, cover);
    Bias_tempa = Phit(tau_b, delta_t)*(bM(tau_a, delta_t) + bT(tau_a, T));
    MSE_tempa  = Phit(tau_b, delta_t)*(VarM(tau_a,delta_t) + VarT(tau_a,T) + ...
                 (bM(tau_a,delta_t)+bT(tau_a,T))^2);
             
    if MissNb > 0
        k = 1:MissNb;  
        tau_k = (k + 1) .* tau_a; 
        
        porb_val = Porb1(tau_b, delta_t);
        phit_val = Phit(tau_b, delta_t);
        porb_k = porb_val .^ k; 
        
        bm_vec = bM(tau_k, delta_t);
        bt_vec = bT(tau_k, T);
        varm_vec = VarM(tau_k, delta_t);
        vart_vec = VarT(tau_k, T);
        
        bias_vec = porb_k .* phit_val .* (bm_vec + bt_vec + k .* tau_a);
        MSE_vec  = porb_k .* phit_val .* (varm_vec + vart_vec + ...
                  (bm_vec + bt_vec).^2 + 2 .* k .* tau_a .* (bm_vec + bt_vec) + (k .* tau_a).^2);
                  
        Bias_tempa = Bias_tempa + sum(bias_vec);
        MSE_tempa  = MSE_tempa  + sum(MSE_vec);
    end
    
    if MissNb ~= 0
        porb_last = Porb1(tau_b, delta_t)^(MissNb+1);
        tau_last  = (MissNb+2)*tau_a;
        bm_last   = bM(tau_last, delta_t);
        bt_last   = bT(tau_last, T);
        
        bias_n1 = porb_last * (bm_last + bt_last + (MissNb+1)*tau_a);
        MSE_n1  = porb_last * (VarM(tau_last,delta_t) + VarT(tau_last,T) + ...
                  (bm_last + bt_last)^2 + ...
                  2*(MissNb+1)*tau_a*(bm_last + bt_last) + ...
                  (MissNb+1)^2*tau_a^2);
                  
        Bias_tempa = Bias_tempa + bias_n1;
        MSE_tempa  = MSE_tempa  + MSE_n1;
    end
    
    Bias_a = Bias_tempa;
    Std_a  = sqrt(max(0, (MSE_tempa - Bias_tempa^2)/N));
    
    %%  tau_b 
    MissNa = missN(tau_a, delta_t, cover);
    Bias_tempb = Phit(tau_a, delta_t)*(bM(tau_b, delta_t) + bT(tau_b, T));
    MSE_tempb  = Phit(tau_a, delta_t)*(VarM(tau_b,delta_t) + VarT(tau_b,T) + ...
                 (bM(tau_b,delta_t)+bT(tau_b,T))^2);
             
    if MissNa > 0
        k = 1:MissNa;  
        tau_k = (k + 1) .* tau_b; 
        
        porb_val = Porb1(tau_a, delta_t);
        phit_val = Phit(tau_a, delta_t);
        porb_k = porb_val .^ k; 
        
        bm_vec = bM(tau_k, delta_t);
        bt_vec = bT(tau_k, T);
        varm_vec = VarM(tau_k, delta_t);
        vart_vec = VarT(tau_k, T);
        
        bias_vec = porb_k .* phit_val .* (bm_vec + bt_vec + k .* tau_b);
        MSE_vec  = porb_k .* phit_val .* (varm_vec + vart_vec + ...
                  (bm_vec + bt_vec).^2 + 2 .* k .* tau_b .* (bm_vec + bt_vec) + (k .* tau_b).^2);
                  
        Bias_tempb = Bias_tempb + sum(bias_vec);
        MSE_tempb  = MSE_tempb  + sum(MSE_vec);
    end
    
    if MissNa ~= 0
        porb_last = Porb1(tau_a, delta_t)^(MissNa+1);
        tau_last  = (MissNa+2)*tau_b;
        bm_last   = bM(tau_last, delta_t);
        bt_last   = bT(tau_last, T);
        
        bias_n1 = porb_last * (bm_last + bt_last + (MissNa+1)*tau_b);
        MSE_n1  = porb_last * (VarM(tau_last,delta_t) + VarT(tau_last,T) + ...
                  (bm_last + bt_last)^2 + ...
                  2*(MissNa+1)*tau_b*(bm_last + bt_last) + ...
                  (MissNa+1)^2*tau_b^2);
                  
        Bias_tempb = Bias_tempb + bias_n1;
        MSE_tempb  = MSE_tempb  + MSE_n1;
    end
    
    Bias_b = Bias_tempb;
    Std_b  = sqrt(max(0, (MSE_tempb - Bias_tempb^2)/N));
    

    Rerr_total   = z_alpha*Std_a/tau_a*100 + z_alpha*Std_b/tau_b*100;
    Rebias_total = abs(Bias_a)/tau_a*100   + abs(Bias_b)/tau_b*100;
end


function y = Porb1(tau, delta_t)
    y = 1 + tau ./ delta_t .* (exp(-delta_t ./ tau) - 1);
end

function y = Phit(tau, delta_t)
    y = tau ./ delta_t .* (1 - exp(-delta_t ./ tau));
end

function y = missN(tau, delta_t, cover)
    y = fix(log(1-cover) ./ log(1 + tau ./ delta_t .* (exp(-delta_t ./ tau)-1))) - 1;
    if y < 0, y = 0; end
end