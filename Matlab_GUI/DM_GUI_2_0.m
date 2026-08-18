function DM_GUI()
% DM MarkovAnalysisGUI  —  Two-State Markov Process Error Analysis (DyeCycling Mode)
% How to run: Type "DM_GUI" directly into the MATLAB command window and press Enter.
%
% NEW: The experiment scale can be specified in TWO ways (choose one):
%   (1) Total number of sampling points N_pts, or
%   (2) Dwell A average count.
% In DM mode the A and B counts are NOT equal; they keep a fixed ratio Na:Nb.
% So the user enters ONE count (interpreted as Dwell A count) and the Dwell B
% count is derived automatically from the fixed ratio.

    %% main window 
    fig = uifigure('Name', 'Two-State DM Markov Process — Uncertainty Analysis', ...
                   'Position', [200 80 500 740], 'Resize', 'off', ...
                   'Color', [0.97 0.97 0.97]);
                   
    uilabel(fig, 'Text', 'DM Two-State Markov Process Uncertainty Analysis', ...
            'Position', [20 690 460 35], 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'FontColor', [0.15 0.25 0.55]);
            
    uipanel(fig, 'Position', [20 680 460 2], 'BackgroundColor', [0.7 0.7 0.8], 'BorderType', 'none');
    
    %% Enter the area label
    lblW = 230; fldX = 270; fldW = 140; rowH = 28; startY = 640; gap = 40;
    
    % tau_a, tau_b
    uilabel(fig, 'Text', '\tau_a (s):', 'Interpreter', 'tex', 'Position', [30 startY lblW rowH], 'FontSize', 12);
    fld_taua = uieditfield(fig, 'numeric', 'Position', [fldX startY+2 fldW 24], 'Value', 1);
    
    uilabel(fig, 'Text', '\tau_b (s):', 'Interpreter', 'tex', 'Position', [30 startY-gap lblW rowH], 'FontSize', 12);
    fld_taub = uieditfield(fig, 'numeric', 'Position', [fldX startY-gap+2 fldW 24], 'Value', 2);
    
    % T_on, T_off
    uilabel(fig, 'Text', 'T_{on} (s):', 'Interpreter', 'tex', 'Position', [30 startY-2*gap lblW rowH], 'FontSize', 12);
    fld_Ton = uieditfield(fig, 'numeric', 'Position', [fldX startY-2*gap+2 fldW 24], 'Value', 4);
    
    uilabel(fig, 'Text', 'T_{off} (s):', 'Interpreter', 'tex', 'Position', [30 startY-3*gap lblW rowH], 'FontSize', 12);
    fld_Toff = uieditfield(fig, 'numeric', 'Position', [fldX startY-3*gap+2 fldW 24], 'Value', 0.5);
    
    % f0
    uilabel(fig, 'Text', 'Sampling frequency f_0 (Hz):', 'Position', [30 startY-4*gap lblW rowH], 'Interpreter', 'tex');
    fld_f0 = uieditfield(fig, 'numeric', 'Position', [fldX startY-4*gap+2 fldW 24], 'Value', 100);

    %% ── Input mode selector ────────────────────────────────────────────────
    uilabel(fig, 'Text', 'Choose ONE input:', ...
        'Position', [30 startY-5.2*gap 440 rowH], 'FontSize', 11, ...
        'FontWeight', 'bold', 'FontColor', [0.15 0.25 0.55]);

    bg = uibuttongroup(fig, ...
        'Position',        [30 startY-6*gap 440 28], ...
        'BorderType',      'none', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'SelectionChangedFcn', @(~,~) updateInputMode());

    rb_Npts = uiradiobutton(bg, 'Text', 'Total sampling points (N)', ...
        'Position', [0 3 210 22], 'FontSize', 11);
    rb_count = uiradiobutton(bg, 'Text', 'Dwell A average count', ...
        'Position', [230 3 210 22], 'FontSize', 11);

    % N_pts field
    lbl_Npts = uilabel(fig, 'Text', 'Total number of sampling points N pts:', ...
        'Position', [30 startY-7*gap lblW rowH]);
    fld_Npts = uieditfield(fig, 'numeric', 'Position', [fldX startY-7*gap+2 fldW 24], 'Value', 1000000);

    % Dwell A count field
    lbl_count = uilabel(fig, 'Text', 'Dwell A average count:', ...
        'Position', [30 startY-7*gap lblW rowH]);
    fld_count = uieditfield(fig, 'numeric', 'Position', [fldX startY-7*gap+2 fldW 24], ...
        'Value', 2180, 'Limits', [1e-9 Inf]);

    % note under the count field
    lbl_note = uilabel(fig, ...
        'Text', 'Note: enter the Dwell A count. Dwell B is derived from the fixed Na:Nb ratio.', ...
        'Position', [30 startY-7.7*gap 450 20], 'FontSize', 10, ...
        'FontColor', [0.45 0.45 0.45], 'FontAngle', 'italic');

    % default: N_pts mode
    rb_Npts.Value = true;
    updateInputMode();

    function updateInputMode()
        if rb_Npts.Value
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
    btn = uibutton(fig, 'Text', '▶  Calculate', 'Position', [165 293 170 38], ...
                   'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.15 0.35 0.75], ...
                   'FontColor', 'white', 'ButtonPushedFcn', @(~,~) runCalc());
    
    %% Result Panel
    resPanel = uipanel(fig, ...
        'Title',           'Analysis Results', ...
        'Position',        [20 15 460 270], ...
        'FontSize',        11, ...
        'FontWeight',      'bold', ...
        'BackgroundColor', [0.94 0.97 1.0], ...
        'BorderType',      'line');
        
    lbl_result = uilabel(resPanel, 'Text', '— Click "Calculate" to view the results —', ...
                         'Position', [15 10 430 230], 'FontSize', 11, 'WordWrap', 'on', ...
                         'Interpreter', 'tex', ...
                         'VerticalAlignment', 'top');
    
    %% Core Computing Logic
    function runCalc()
        btn.Enable = 'off'; 
        lbl_result.Text = 'Calculating...'; 
        lbl_result.FontColor = [0.4 0.4 0.4];
        drawnow;
        
        try
            tau_a = fld_taua.Value;  tau_b = fld_taub.Value;
            T_on  = fld_Ton.Value;   T_off = fld_Toff.Value;
            f0    = fld_f0.Value;

            % Determine N_pts from whichever input mode is active
            if rb_Npts.Value
                N_pts = fld_Npts.Value;
            else
                countA = fld_count.Value;
                N_pts  = Npts_from_countA(countA, tau_a, tau_b, T_on, T_off, f0);
            end

            [Rebias_total, Rerr_total, Bias_Total_a, Bias_Total_b, ...
             Std_Total_a, Std_Total_b, M, Na, Nb] = ...
                computeDM(tau_a, tau_b, T_on, T_off, f0, N_pts);

            % Derived counterpart line for cross-checking
            if rb_Npts.Value
                derived_str = sprintf('Derived Dwell A count : %.1f   |   N pts : %.0f', Na*M, N_pts);
            else
                derived_str = sprintf('Derived total sampling points N : %.0f', N_pts);
            end

            resStr = { ...
                sprintf('Total Relative Bias: %.4f %%', Rebias_total), ...
                sprintf('Total Relative Variation: %.4f %%', Rerr_total), ...
                '-----------------------------------------------------------------', ...
                sprintf('\\tau_a Bias : %.4e ', Bias_Total_a), ...
                sprintf('\\tau_a Variation : %.4e ', Std_Total_a), ...
                sprintf('\\tau_b Bias : %.4e ', Bias_Total_b), ...
                sprintf('\\tau_b Variation : %.4e ', Std_Total_b), ...
                '-----------------------------------------------------------------', ...
                sprintf('Dwell A average count    : %.1f',    Na * M), ...
                sprintf('Dwell B average count    : %.1f',    Nb * M), ...
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
end

%%  Invert Dwell A count - N_pts
% In DM mode Na, Nb do not depend on T, so the inversion is direct and linear:
%   M = countA / Na ;  T = M*(T_on+T_off) ;  N_pts = T*f0
function N_pts = Npts_from_countA(countA, tau_a, tau_b, T_on, T_off, f0)
    [Na, ~] = computeNaNb(tau_a, tau_b, T_on, T_off);
    M = countA / Na;
    T = M * (T_on + T_off);
    N_pts = T * f0;
end

% Na, Nb evaluated with E = 0 (large-T limit), matching the values used in
% computeDM for any realistic experiment length.
function [Na, Nb] = computeNaNb(tau_a, tau_b, T_on, T_off)
    A = tau_a/(tau_a+tau_b); B = tau_b/(tau_a+tau_b);
    E = 0;                                     % large-T limit
    P1 = A*(A+B*E); P2 = A*(B-B*E); P3 = B*(A-A*E); P4 = B*(B+A*E);
    TA = (2*P1+1.5*P2+1.5*P3+P4)*(tau_a+tau_b);
    TB = (2*P4+1.5*P2+1.5*P3+P1)*(tau_a+tau_b);
    Na = (TA+T_on)/(tau_a+tau_b)*exp(-TA/T_on);
    Nb = (TB+T_on)/(tau_a+tau_b)*exp(-TB/T_on);
end

%% Core computation (unchanged from original)
function [Rebias_total, Rerr_total, Bias_Total_a, Bias_Total_b, Std_Total_a, Std_Total_b, M, Na, Nb] = ...
         computeDM(tau_a, tau_b, T_on, T_off, f0, N_pts)
         
    alpha = 0.05; z_alpha = -norminv(alpha/2);
    delta_t = 1/f0;
    T = N_pts / f0; 
    cover = 0.9;
    
    M = max(1, T/(T_on + T_off));
    Step = 500;
    
  
    A = tau_a/(tau_a+tau_b); B = tau_b/(tau_a+tau_b);
    E = exp(-T*(tau_a+tau_b)/(tau_a*tau_b));
    P1 = A*(A+B*E); P2 = A*(B-B*E); P3 = B*(A-A*E); P4 = B*(B+A*E);
    TA = (2*P1+1.5*P2+1.5*P3+P4)*(tau_a+tau_b);
    TB = (2*P4+1.5*P2+1.5*P3+P1)*(tau_a+tau_b);
    
 
    xa_edges = (10.^(0.01.*(0:Step)))*1e-2;
    xb_edges = (10.^(0.01.*(0:Step)))*1e-2;
    xa = xa_edges(1:Step); dxa = diff(xa_edges);
    xb = xb_edges(1:Step); dxb = diff(xb_edges);
    
    Na = (TA+T_on)/(tau_a+tau_b)*exp(-TA/T_on);
    Nb = (TB+T_on)/(tau_a+tau_b)*exp(-TB/T_on);
    Ifactora = exp(-xa/T_on)/T_on .* dxa;
    Ifactorb = exp(-xb/T_on)/T_on .* dxb;
    

    MissNb = missN(tau_b, delta_t, cover);
    
  
    Bias_tempa = Phit(tau_b, delta_t) .* (bM(tau_a, delta_t) + bT(tau_a, xa));
    MSE_tempa  = Phit(tau_b, delta_t) .* (VarM(tau_a, delta_t) + VarT(tau_a, xa) + (bM(tau_a, delta_t) + bT(tau_a, xa)).^2);
    
    if MissNb > 0
   
        k_vec = (1:MissNb)'; 
        tau_k = (k_vec + 1) * tau_a;
        [TAU_mat, XA_mat] = ndgrid(tau_k, xa);
        
        PORB_mat = repmat(Porb1(tau_b, delta_t).^k_vec, 1, Step);
        BM_mat   = repmat(bM(tau_k, delta_t), 1, Step);
        K_mat    = repmat(k_vec, 1, Step);
        
        BT_mat = bT(TAU_mat, XA_mat);
        VT_mat = VarT(TAU_mat, XA_mat);
        phit_val = Phit(tau_b, delta_t);
        
  
        bias_n_mat = PORB_mat .* phit_val .* (BM_mat + BT_mat + K_mat * tau_a);
        MSE_n_mat  = PORB_mat .* phit_val .* (VarM(TAU_mat, delta_t) + VT_mat + ...
                     (BM_mat + BT_mat).^2 + 2 .* K_mat .* tau_a .* (BM_mat + BT_mat) + (K_mat .* tau_a).^2);
                     
        Bias_tempa = Bias_tempa + sum(bias_n_mat, 1);
        MSE_tempa  = MSE_tempa  + sum(MSE_n_mat, 1);
    end
    
    if MissNb ~= 0
        tau_last = (MissNb+2)*tau_a;
        PORB_last = Porb1(tau_b, delta_t)^(MissNb+1);
        BM_last = bM(tau_last, delta_t);
        
        bias_n1 = PORB_last .* (BM_last + bT(tau_last, xa) + (MissNb+1)*tau_a);
        MSE_n1  = PORB_last .* (VarM(tau_last, delta_t) + VarT(tau_last, xa) + ...
                  (BM_last + bT(tau_last, xa)).^2 + 2*(MissNb+1)*tau_a.*(BM_last + bT(tau_last, xa)) + ((MissNb+1)*tau_a)^2);
                  
        Bias_tempa = Bias_tempa + bias_n1;
        MSE_tempa  = MSE_tempa  + MSE_n1;
    end
    

    temp1_a = (MSE_tempa - Bias_tempa.^2) ./ Na;
    IVTa = sum(temp1_a .* Ifactora);
    IBTa = sum(Bias_tempa .* Ifactora);


    MissNa = missN(tau_a, delta_t, cover);
    
    Bias_tempb = Phit(tau_a, delta_t) .* (bM(tau_b, delta_t) + bT(tau_b, xb));
    MSE_tempb  = Phit(tau_a, delta_t) .* (VarM(tau_b, delta_t) + VarT(tau_b, xb) + (bM(tau_b, delta_t) + bT(tau_b, xb)).^2);
    
    if MissNa > 0
        k_vec = (1:MissNa)'; 
        tau_k = (k_vec + 1) * tau_b;
        [TAU_mat, XB_mat] = ndgrid(tau_k, xb);
        
        PORB_mat = repmat(Porb1(tau_a, delta_t).^k_vec, 1, Step);
        BM_mat   = repmat(bM(tau_k, delta_t), 1, Step);
        K_mat    = repmat(k_vec, 1, Step);
        
        BT_mat = bT(TAU_mat, XB_mat);
        VT_mat = VarT(TAU_mat, XB_mat);
        phit_val = Phit(tau_a, delta_t);
        
        bias_n_mat = PORB_mat .* phit_val .* (BM_mat + BT_mat + K_mat * tau_b);
        MSE_n_mat  = PORB_mat .* phit_val .* (VarM(TAU_mat, delta_t) + VT_mat + ...
                     (BM_mat + BT_mat).^2 + 2 .* K_mat .* tau_b .* (BM_mat + BT_mat) + (K_mat .* tau_b).^2);
                     
        Bias_tempb = Bias_tempb + sum(bias_n_mat, 1);
        MSE_tempb  = MSE_tempb  + sum(MSE_n_mat, 1);
    end
    
    if MissNa ~= 0
        tau_last = (MissNa+2)*tau_b;
        PORB_last = Porb1(tau_a, delta_t)^(MissNa+1);
        BM_last = bM(tau_last, delta_t);
        
        bias_n1 = PORB_last .* (BM_last + bT(tau_last, xb) + (MissNa+1)*tau_b);
        MSE_n1  = PORB_last .* (VarM(tau_last, delta_t) + VarT(tau_last, xb) + ...
                  (BM_last + bT(tau_last, xb)).^2 + 2*(MissNa+1)*tau_b.*(BM_last + bT(tau_last, xb)) + ((MissNa+1)*tau_b)^2);
                  
        Bias_tempb = Bias_tempb + bias_n1;
        MSE_tempb  = MSE_tempb  + MSE_n1;
    end
    
    temp1_b = (MSE_tempb - Bias_tempb.^2) ./ Nb;
    IVTb = sum(temp1_b .* Ifactorb);
    IBTb = sum(Bias_tempb .* Ifactorb);

    % I0a = sum(Ifactora);   % = I_{0,DM} 
    % I0b = sum(Ifactorb);

    Var_Total_a = IVTa/M;  
    Bias_Total_a = IBTa;
    Var_Total_b = IVTb/M;  
    Bias_Total_b = IBTb;
    Std_Total_a = sqrt(max(0, Var_Total_a));
    Std_Total_b = sqrt(max(0, Var_Total_b));
    
    Rerr_total   = z_alpha*Std_Total_a/tau_a*100 + z_alpha*Std_Total_b/tau_b*100;
    Rebias_total = abs(Bias_Total_a)/tau_a*100 + abs(Bias_Total_b)/tau_b*100;
end


function y = Porb1(tau, delta_t)
   y = 1 + tau./delta_t.*(exp(-delta_t./tau)-1);
end

function y = Phit(tau, delta_t)
    y = tau./delta_t.*(1-exp(-delta_t./tau));
end

function y = missN(tau, delta_t, cover)
    y = fix(log(1-cover)./log(1+tau./delta_t.*(exp(-delta_t./tau)-1))) - 1;
    if y < 0, y = 0; end
end

function y = bM(tau, delta_t)
   A = exp(-delta_t./tau);
   y = delta_t.*(A+1)./(1-A) - 2.*tau;
end

function y = bT(tau, T)

   if isscalar(tau) && ~isscalar(T), tau = repmat(tau, size(T)); end
   if ~isscalar(tau) && isscalar(T), T = repmat(T, size(tau)); end
   
   y = 2.*T.^2./tau.*exp(T./tau).*(expint(2.*T./tau)-expint(T./tau)) + (T-tau).*(2-exp(-T./tau))+tau;
   
   bad_idx = isnan(y) | isinf(y);
   if any(bad_idx(:))
        t_bad = tau(bad_idx); T_bad = T(bad_idx);
        y(bad_idx) = exp(-T_bad./t_bad).*(t_bad./2+t_bad.^2./(2.*T_bad)) - 4.*t_bad.^2./T_bad;
   end
   
   I1a = 2.*T./tau.*exp(T./tau).*(igamma(0, T./tau)-igamma(0, 2.*T./tau)) + exp(-T./tau) - 1;
   bad_i = isnan(I1a) | isinf(I1a);
   if any(bad_i(:)), I1a(bad_i) = 1; end
   
   y = y ./ I1a-tau;
end

function y = VarM(tau, delta_t)
    A = exp(-delta_t./tau);
    y = delta_t.^2.*A ./ (6.*tau./delta_t.*(1-A));
end

function y = VarT(tau, T)
   if isscalar(tau) && ~isscalar(T), tau = repmat(tau, size(T)); end
   if ~isscalar(tau) && isscalar(T), T = repmat(T, size(tau)); end

   y = (-2.*T.^2./tau.*exp(T./tau).*(expint(2.*T./tau)-expint(T./tau)).*(T+2.*tau) + (exp(-T./tau)-1).*(T.^2+2.*T.*tau)-T.^2+tau.^2) ...
       - (2.*T.^2./tau.*exp(T./tau).*(expint(2.*T./tau)-expint(T./tau)) + (T-tau).*(2-exp(-T./tau))).^2;
            
   bad_y = isnan(y) | isinf(y);
   if any(bad_y(:)), y(bad_y) = tau(bad_y).^2; end
   
   I1a = 2.*T./tau.*exp(T./tau).*(igamma(0, T./tau)-igamma(0, 2.*T./tau)) + exp(-T./tau) - 1;
   bad_i = isnan(I1a) | isinf(I1a);
   if any(bad_i(:)), I1a(bad_i) = 1; end
   
   y = y ./ I1a;    
end