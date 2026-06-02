% =========================================================================
% TPTL 稳态分布参数模型 (异构变径并联环路专用版 - 彻底全域平滑高精度版)
% =========================================================================
clear; clc; close all;

%% 1. 几何参数与边界条件
geo.L_down_v = 3.0;  geo.N_down_v = 30;  geo.dz_down_v = geo.L_down_v / geo.N_down_v;
geo.D_i_down = 0.05; geo.A_c_down = pi/4 * geo.D_i_down^2;

geo.L_down_h = 3.0;  geo.N_down_h = 30;  geo.dz_down_h = geo.L_down_h / geo.N_down_h;

geo.N_tubes_evap = 3;
geo.L_evap = 2.83;   geo.N_evap = 50;    geo.dz_evap = geo.L_evap / geo.N_evap;
geo.D_i_evap = 0.02; geo.D_o_evap = 0.025; geo.A_c_evap = pi/4 * geo.D_i_evap^2;

geo.L_riser_v = 1.2; geo.N_riser_v = 20; geo.dz_riser_v = geo.L_riser_v / geo.N_riser_v;
geo.D_i_riser = 0.05; geo.A_c_riser = pi/4 * geo.D_i_riser^2;

geo.N_tubes_cond = 3;
geo.L_cond = 3.0;    geo.N_cond = 50;    geo.dz_cond = geo.L_cond / geo.N_cond;
geo.D_i_cond = 0.02; geo.D_o_cond = 0.025; geo.A_c_cond = pi/4 * geo.D_i_cond^2;

geo.k_steel = 22.0;  
geo.K_bend  = 1.2;   

geo.T_wall_evap = 273.15 + 250; 
geo.T_wall_cond = 273.15 + 20; 
M_charge = 10; % [kg] 系统初始充注量

%% 2. 求解器初始猜测值与边界设置
X0 = [2, 1.5, 0.1];        
lb = [0.000611657, 0.02, 1e-4];
P_max = safe_IAPWS('psat_T', geo.T_wall_evap); 
ub = [P_max, geo.L_down_v, 5];

%% 3. 配置并调用 MultiStart 全局求解器
options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective', ...
    'FunctionTolerance', 1e-8, 'StepTolerance', 1e-10, 'MaxIterations', 500); 
problem = createOptimProblem('lsqnonlin', 'x0', X0, ...
    'objective', @(X) TPTL_residuals(X, geo, M_charge, X0), 'lb', lb, 'ub', ub, 'options', options);
ms = MultiStart('Display', 'iter', 'UseParallel', false);
fprintf('开始调用 MultiStart 全局搜索 (启动全域平滑物性内核)...\n');
tic;
[X_sol, resnorm, ~, ~, solutions] = run(ms, problem, 10); 
toc;

F_final = TPTL_residuals(X_sol, geo, M_charge, X0);
fprintf('\n=== 全局求解完成 ===\n');
fprintf('最优解: P0 = %.4f MPa, H0 = %.4f m, 总流量 W = %.4f kg/s\n', X_sol(1), X_sol(2), X_sol(3));
fprintf('最小残差平方和 (resnorm) = %e\n', resnorm);
fprintf('压力方程相对误差 (err_P) = %+.4e\n', F_final(1));
fprintf('比焓方程相对误差 (err_h) = %+.4e\n', F_final(2));
fprintf('质量方程相对误差 (err_M) = %+.4e\n', F_final(3));

%% 4. 后处理调用
fprintf('\n==================================================\n');
fprintf('>>> 正在全自动调用后处理函数 TPTL_postprocess ...\n');
TPTL_postprocess_2026_5_30_all_smooth(X_sol(1), X_sol(2), X_sol(3));

%% ===================== 局部函数区 (平滑版) =====================
function F = TPTL_residuals(X, geo, M_charge, X0)
    P0 = X(1); H0 = X(2); W = X(3); 
    try
        h0 = safe_IAPWS('hL_p', P0);
        P_curr = P0; h_curr = h0; 
        Total_Mass = 0;
        R_wall_evap = (geo.D_i_evap / (2 * geo.k_steel)) * log(geo.D_o_evap / geo.D_i_evap);
        R_wall_cond = (geo.D_i_cond / (2 * geo.k_steel)) * log(geo.D_o_cond / geo.D_i_cond);
        
        G_down  = W / geo.A_c_down; G_evap  = (W / geo.N_tubes_evap) / geo.A_c_evap;
        G_riser = W / geo.A_c_riser; G_cond  = (W / geo.N_tubes_cond) / geo.A_c_cond;
        
        % (1) 下降段 - 竖直 
        L_unfilled = max(0, geo.L_down_v - H0); curr_z = 0;
        for i = 1:geo.N_down_v
            z_start = curr_z; curr_z = curr_z + geo.dz_down_v; z_end = curr_z;
            if L_unfilled >= z_end
                rho = 1 / safe_IAPWS('vV_p', P_curr); 
                Total_Mass = Total_Mass + rho * geo.A_c_down * geo.dz_down_v;
            elseif L_unfilled <= z_start
                [rho, mu, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
                Re = G_down * geo.D_i_down / mu; f = calc_fanning_friction(Re);      
                dP_f = f * (geo.dz_down_v/geo.D_i_down) * (G_down^2)/(2*rho);
                dP_g = rho * 9.81 * geo.dz_down_v; 
                P_curr = P_curr + ((dP_g - dP_f) / 1e6); 
                Total_Mass = Total_Mass + rho * geo.A_c_down * geo.dz_down_v;
            else
                L_gas = L_unfilled - z_start; L_liq = z_end - L_unfilled;
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mass_g = rho_g * geo.A_c_down * L_gas;
                [rho_l, mu_l, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
                Re = G_down * geo.D_i_down / mu_l; f_l = calc_fanning_friction(Re);
                dP_f_l = f_l * (L_liq/geo.D_i_down) * (G_down^2)/(2*rho_l); dP_g_l = rho_l * 9.81 * L_liq;
                Total_Mass = Total_Mass + mass_g + rho_l * geo.A_c_down * L_liq; 
                P_curr = P_curr + ((dP_g_l - dP_f_l) / 1e6);
            end
        end
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_down, geo.K_bend) / 1e6;
        
        % (2) 下降段 - 水平 
        for i = 1:geo.N_down_h
            [~, mu_s, ~, x_th, ~, rho_mix, ~] = get_fluid_state(P_curr, h_curr);
            [rho_l, mu_l, rho_g, mu_g, ~, ~, ~, ~] = get_safe_props(P_curr, h_curr, x_th);
            x_safe_tp = max(min(x_th, 0.999), 1e-4);
            w_liq = 0.5 * (1 - tanh((x_th - 0.001) / 0.001)); w_gas = 0.5 * (1 + tanh((x_th - 0.999) / 0.001)); w_tp = max(1 - w_liq - w_gas, 0);
            
            Re_sp = G_down * geo.D_i_down / mu_s; f_sp = calc_fanning_friction(Re_sp); 
            dP_f_sp = f_sp * (geo.dz_down_h/geo.D_i_down) * (G_down^2)/(2*rho_mix);
            dP_f_tp = calc_two_phase_friction_dp(x_safe_tp, G_down, geo.D_i_down, rho_l, mu_l, rho_g, mu_g, geo.dz_down_h);
            
            P_curr = P_curr - ((w_liq * dP_f_sp + w_tp * dP_f_tp + w_gas * dP_f_sp) / 1e6); Total_Mass = Total_Mass + rho_mix * geo.A_c_down * geo.dz_down_h;
        end
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_down, geo.K_bend) / 1e6;
        
        % (3) 蒸发段 - 竖直 
        for i = 1:geo.N_evap
            [~, ~, T_fluid, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            [rho_l, mu_l, rho_g, mu_g, Cp_fluid, k_fluid, ~, ~] = get_safe_props(P_curr, h_curr, x_th);
            x_safe_tp = max(min(x_th, 0.999), 1e-4);
            w_liq = 0.5 * (1 - tanh((x_th - 0.001) / 0.001)); w_gas = 0.5 * (1 + tanh((x_th - 0.999) / 0.001)); w_tp = max(1 - w_liq - w_gas, 0);
            
            Re_sp = G_evap * geo.D_i_evap / mu_mix;
            Pr_sp = (Cp_fluid * 1000) * mu_mix / k_fluid;
            h_wall = safe_IAPWS('h_pT', P_curr, geo.T_wall_evap); mu_wall = safe_IAPWS('mu_ph', P_curr, h_wall); 
            h_htc_sp = calc_single_phase_htc(Re_sp, Pr_sp, k_fluid, geo.D_i_evap, max(i * geo.dz_evap, 1e-4), mu_mix, mu_wall, true);
            q_node_sp = (1 / (1 / h_htc_sp + R_wall_evap)) * (geo.T_wall_evap - T_fluid);
            
            % q_guess = 5000; tol_q = 1e-3;
            % for iter = 1:50 
            %     q_new = calc_boiling_htc(P_curr, q_guess, x_safe_tp) * ((geo.T_wall_evap - q_guess * R_wall_evap) - T_fluid); 
            %     if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end; q_guess = 0.8 * q_guess + 0.2 * q_new; 
            % end
            h_htc=calc_boiling_htc(P_curr, 9000000, x_safe_tp);
          
            q_new = (1 / (1 / h_htc + R_wall_evap)) * (geo.T_wall_evap - T_fluid);
            dQ_total = geo.N_tubes_evap * ((w_liq * q_node_sp + w_tp * q_new + w_gas * q_node_sp) * pi * geo.D_i_evap * geo.dz_evap);
     
            dh = (dQ_total / W) / 1000; dP_g = rho_mix * 9.81 * geo.dz_evap;
            
            f_sp = calc_fanning_friction(Re_sp); dP_f_sp = f_sp * (geo.dz_evap/geo.D_i_evap) * (G_evap^2)/(2*rho_mix);
            dP_f_tp = calc_two_phase_friction_dp(x_safe_tp, G_evap, geo.D_i_evap, rho_l, mu_l, rho_g, mu_g, geo.dz_evap);
            
            P_curr = P_curr - ((dP_g + w_liq * dP_f_sp + w_tp * dP_f_tp + w_gas * dP_f_sp) / 1e6); h_curr = h_curr + dh; 
            Total_Mass = Total_Mass + geo.N_tubes_evap * rho_mix * geo.A_c_evap * geo.dz_evap;
        end
        
        % (4) 上升段 - 竖直 
        for i = 1:geo.N_riser_v
            [~, ~, ~, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            [rho_l, mu_l, rho_g, mu_g, ~, ~, ~, ~] = get_safe_props(P_curr, h_curr, x_th);
            x_safe_tp = max(min(x_th, 0.999), 1e-4);
            w_liq = 0.5 * (1 - tanh((x_th - 0.001) / 0.001)); w_gas = 0.5 * (1 + tanh((x_th - 0.999) / 0.001)); w_tp = max(1 - w_liq - w_gas, 0);
            
            Re_sp = G_riser * geo.D_i_riser / mu_mix; f_sp = calc_fanning_friction(Re_sp); 
            dP_f_sp = f_sp * (geo.dz_riser_v/geo.D_i_riser) * (G_riser^2)/(2*rho_mix);
            dP_f_tp = calc_two_phase_friction_dp(x_safe_tp, G_riser, geo.D_i_riser, rho_l, mu_l, rho_g, mu_g, geo.dz_riser_v);
            
            dP_g = rho_mix * 9.81 * geo.dz_riser_v;
            P_curr = P_curr - ((dP_g + w_liq * dP_f_sp + w_tp * dP_f_tp + w_gas * dP_f_sp) / 1e6); 
            Total_Mass = Total_Mass + rho_mix * geo.A_c_riser * geo.dz_riser_v;
        end
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_riser, geo.K_bend) / 1e6;
        
        % (5) 冷凝段 - 水平 
        for i = 1:geo.N_cond
            [~, ~, T_fluid, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            [rho_l, mu_l, rho_g, mu_g, Cp_fluid, k_fluid, k_l, h_fg] = get_safe_props(P_curr, h_curr, x_th);
            x_safe_tp = max(min(x_th, 0.999), 1e-4);
            w_liq = 0.5 * (1 - tanh((x_th - 0.001) / 0.001)); w_gas = 0.5 * (1 + tanh((x_th - 0.999) / 0.001)); w_tp = max(1 - w_liq - w_gas, 0);
            
            Re_sp = G_cond * geo.D_i_cond / mu_mix; Pr_sp = (Cp_fluid * 1000) * mu_mix / k_fluid; 
            h_wall = safe_IAPWS('h_pT', P_curr, max(T_fluid - 10, 273.16)); mu_wall = safe_IAPWS('mu_ph', P_curr, h_wall); 
            z_curr = max(i * geo.dz_cond, 1e-4); 
            h_htc_sp = calc_single_phase_htc(Re_sp, Pr_sp, k_fluid, geo.D_i_cond, z_curr, mu_mix, mu_wall, false);
            q_node_sp = (1 / (1 / h_htc_sp + R_wall_cond)) * max(T_fluid - geo.T_wall_cond, 0.1);
            
            q_guess = 5000; tol_q = 1e-3;
            for iter = 1:50
                dT_film = max(T_fluid - (geo.T_wall_cond + q_guess * R_wall_cond), 0.1); 
                q_new = 1.13 * ( (k_l^3 * rho_l^2 * h_fg * 9.81) / (mu_l * z_curr * dT_film) )^0.25 * dT_film;
                if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end; q_guess = 0.8 * q_guess + 0.2 * q_new; 
            end
            
            dQ_total = geo.N_tubes_cond * ((w_liq * q_node_sp + w_tp * q_new + w_gas * q_node_sp) * pi * geo.D_i_cond * geo.dz_cond); 
            dh = -(dQ_total / W) / 1000;
            
            f_sp = calc_fanning_friction(Re_sp); dP_f_sp = f_sp * (geo.dz_cond/geo.D_i_cond) * (G_cond^2)/(2*rho_mix);
            dP_f_tp = calc_two_phase_friction_dp(x_safe_tp, G_cond, geo.D_i_cond, rho_l, mu_l, rho_g, mu_g, geo.dz_cond);
            
            P_curr = P_curr - ((w_liq * dP_f_sp + w_tp * dP_f_tp + w_gas * dP_f_sp) / 1e6); h_curr = h_curr + dh; 
            Total_Mass = Total_Mass + geo.N_tubes_cond * rho_mix * geo.A_c_cond * geo.dz_cond;
        end
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_down, geo.K_bend) / 1e6;
        
        F = [(P_curr - P0) / P0, (h_curr - h0) / h0, (Total_Mass - M_charge) / M_charge];
    catch ME
        F = 1e5 * [(X(1)-X0(1))^2 + 1, (X(2)-X0(2))^2 + 1, (X(3)-X0(3))^2 + 1]; 
    end
end

% ======= 重点物理属性与平滑内核区 =======
function [rho_l, mu_l, rho_g, mu_g, cp_sp, k_sp, k_l, h_fg] = get_safe_props(P, h, x_th)
    % 该函数专门在底层阻断任何将饱和两相直接扔进物性库引发退出的风险
    Tsat = safe_IAPWS('Tsat_p', P);
    T_liq_safe = max(Tsat - 0.1, 273.16); 
    T_gas_safe = Tsat + 0.1;
    
    hL_safe = safe_IAPWS('h_pT', P, T_liq_safe);
    hV_safe = safe_IAPWS('h_pT', P, T_gas_safe);
    
    rho_l = 1 / safe_IAPWS('vL_p', P); rho_g = 1 / safe_IAPWS('vV_p', P);
    mu_l = safe_IAPWS('mu_ph', P, hL_safe); mu_g = safe_IAPWS('mu_ph', P, hV_safe);
    k_l = safe_IAPWS('k_ph', P, hL_safe); h_fg = (safe_IAPWS('hV_p', P) - safe_IAPWS('hL_p', P)) * 1000;
    
    if x_th <= 0
        cp_sp = safe_IAPWS('cp_ph', P, h); k_sp = safe_IAPWS('k_ph', P, h);
    elseif x_th >= 1
        cp_sp = safe_IAPWS('cp_ph', P, h); k_sp = safe_IAPWS('k_ph', P, h);
    else
        % 当处于两相区时，强行利用亚冷液或过热气状态返回单相参数虚假占位
        if x_th < 0.5
            cp_sp = safe_IAPWS('cp_ph', P, hL_safe); k_sp = safe_IAPWS('k_ph', P, hL_safe);
        else
            cp_sp = safe_IAPWS('cp_ph', P, hV_safe); k_sp = safe_IAPWS('k_ph', P, hV_safe);
        end
    end
end

function [rho_single, mu_single, T_fluid, x_th, alpha, rho_mix, mu_mix] = get_fluid_state(P, h)
    hL = safe_IAPWS('hL_p', P); hV = safe_IAPWS('hV_p', P); Tsat = safe_IAPWS('Tsat_p', P); x_th = (h - hL) / (hV - hL);
    if x_th <= 0 
        x_th = 0; T_fluid = safe_IAPWS('T_ph', P, h); rho_single = 1 / safe_IAPWS('v_ph', P, h); mu_single = safe_IAPWS('mu_ph', P, h); alpha = 0; rho_mix = rho_single; mu_mix = mu_single;
    elseif x_th >= 1
        x_th = 1; T_fluid = safe_IAPWS('T_ph', P, h); rho_single = 1 / safe_IAPWS('v_ph', P, h); mu_single = safe_IAPWS('mu_ph', P, h); alpha = 1; rho_mix = rho_single; mu_mix = mu_single;
    else
        T_fluid = Tsat; vL = safe_IAPWS('vL_p', P); vV = safe_IAPWS('vV_p', P); rho_l = 1 / vL; rho_g = 1 / vV;
        T_liq_safe = max(Tsat - 0.1, 273.16); T_gas_safe = Tsat + 0.1;
        mu_l = safe_IAPWS('mu_ph', P, safe_IAPWS('h_pT', P, T_liq_safe)); mu_g = safe_IAPWS('mu_ph', P, safe_IAPWS('h_pT', P, T_gas_safe));
        alpha = max(0, min(1, 1 / (1 + ((1-x_th)/x_th) * (rho_g/rho_l)^0.89 * (mu_l/mu_g)^0.18)));
        rho_mix = alpha * rho_g + (1 - alpha) * rho_l; mu_mix = 1 / ( (x_th/mu_g) + ((1-x_th)/mu_l) ); rho_single = rho_l; mu_single = mu_l;
    end
end

function dP_local = calc_local_dp(P, x, G, K_bend)
    rho_l = 1 / safe_IAPWS('vL_p', P); rho_g = 1 / safe_IAPWS('vV_p', P); 
    dP_l_only = K_bend * (G^2) / (2 * rho_l); dP_g_only = K_bend * (G^2) / (2 * rho_g);
    x_safe = max(min(x, 0.999), 1e-4); Xb = (rho_g/rho_l)^0.5 * (1-x_safe)/x_safe; 
    phi_c_sq = 1 + ((1 + (3.2-1)*(rho_g)^0.5) * ((rho_g/rho_l)^0.5 + (rho_g/rho_l)^-0.5))/Xb + 1/Xb^2; 
    w_liq = 0.5 * (1 - tanh((x - 0.005) / 0.002)); w_gas = 0.5 * (1 + tanh((x - 0.995) / 0.002)); w_tp = max(1 - w_liq - w_gas, 0); 
    dP_local = w_liq * dP_l_only + w_tp * (dP_l_only *(1-x_safe)^2* phi_c_sq) + w_gas * dP_g_only;
end

function f_fan = calc_fanning_friction(Re)
    Re = max(Re, 1e-4); weight = 0.5 * (1 + tanh((Re - 2000) / 100)); 
    f_fan = (1 - weight) * (64 / Re) + weight * (0.184 / Re^0.2);
end

function dP_f = calc_two_phase_friction_dp(x, G, D, rho_l, mu_l, rho_g, mu_g, dz)
    G_l = G * (1 - x); G_g = G * x; Re_l = max(G_l * D / mu_l, 1e-4); Re_g = max(G_g * D / mu_g, 1e-4);
    f_l = calc_fanning_friction(Re_l); X = max(((1 - x) / x) * sqrt((f_l / calc_fanning_friction(Re_g)) * (rho_g / rho_l)), 1e-6);
    w_l = 0.5 * (1 + tanh((Re_l - 2000) / 100)); w_g = 0.5 * (1 + tanh((Re_g - 2000) / 100));
    C = (1 - w_l) * (1 - w_g) * 5 + (1 - w_l) * w_g * 12 + w_l * (1 - w_g) * 10 + w_l * w_g * 20;
    dP_f = (1 + C / X + 1 / X^2) * (f_l * (G_l^2 / rho_l) / (2 * D)) * dz;
end

function h_htc = calc_single_phase_htc(Re, Pr, k, D, z, mu_fluid, mu_wall, is_heating)
    Re = max(Re, 10); Pr = max(Pr, 0.1); z = max(z, 1e-4); weight = 0.5 * (1 + tanh((Re - 2200) / 100));
    n = 0.4; if ~is_heating, n = 0.3; end
    h_htc = (1 - weight) * (1.24 * (k * (1 + 5 * exp(-z / (10 * D))) / D) * (Re * Pr * (D / z))^(1/3) * (mu_fluid / mu_wall)) ...
          + weight * (0.027 * (k * (1 + 7 * (D / z)) / D) * Re^0.8 * Pr^n * (mu_fluid / mu_wall)^0.14);
end

function val = safe_IAPWS(prop, varargin)
    val = IAPWS_IF97(prop, varargin{:}); if any(isnan(val)), error('TPTL:NaN_Error', 'IAPWS_IF97 out of bounds.'); end
end

function h_htc = calc_boiling_htc(P_MPa, q, x)
    h_htc = 540 * max(P_MPa / 22.064, 1e-4)^0.394 *(1-max(0, min(x, 0.99)))^-0.65* 18.015^(-0.5) * max(q, 10)^0.54;
end
