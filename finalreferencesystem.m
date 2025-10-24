% ---------- DATA ----------
filename = 'Measurements Intex Solarmat 2024-07-03.xlsx';
data = readtable(filename,'Sheet',1,'Range','A6:C380');
t_meas   = data{:,1};   Tin_meas = data{:,2};   Tout_meas= data{:,3};

% ---------- GIVEN ----------
Vdot = 1.48/60*1e-3;        % m^3/s
T_env0 = 20 + 273.15;       % K  (far-field room/ground)
T_ini  = 20.4 + 273.15;     % K
G     = 1000;               % W/m^2
eta0  = 0.25;               % -
A_mat = 1.0;                % m^2

% Tank / material
V_wt=2e-3; H_t=0.20; D_t=0.11; d_wall=0.003; d_lid=0.005; d_1=0.05;
k_1=0.04; k_pvc=0.19; rho_w=0.997e3; cp=4184;

% Mat sheet + ground
d_mat=0.0015; k_mat=0.19;
h_conv_air = 5.2;            % W/m^2K (top convection)
h_ground   = 7;           % W/m^2K (bottom contact to floor)
eps_mat    = 0.85;         % [-] mat emissivity (top surface)

% Pipes
L_p=3.0; Di_p=0.010; t_p=0.0015; t_ins_p=0.010; k_pipe=0.19; k_ins=0.04; eps_pipe=0.85;

% Constants
sigmaSB = 5.670374419e-8;

% ---------- DERIVED ----------
mdot = rho_w*Vdot;  m_wt = rho_w*V_wt;

% Tank geometry / base resistances
r1=D_t/2; r2=r1+d_wall; r3=r2+d_1;
A_floor = pi*r1^2; A_wall_outer = 2*pi*r3*H_t; A_lid = A_floor;
R_wall_pvc = log(r2/r1)/(2*pi*k_pvc*H_t);
R_wall_ins = log(r3/r2)/(2*pi*k_1*H_t);
R_cond_wall = R_wall_pvc + R_wall_ins;                  % K/W
R_cond_floor = d_wall/(k_pvc*A_floor) + d_1/(k_1*A_floor);
R_cond_lid   = d_lid /(k_pvc*A_lid);

% Pipes to AIR
ri_p=Di_p/2; ro_p=ri_p+t_p; ro_i=ro_p+t_ins_p; A_pipe_ext=2*pi*ro_i*L_p;
R_pipe_pvc = log(ro_p/ri_p)/(2*pi*k_pipe*L_p);
R_pipe_ins = log(ro_i/ro_p)/(2*pi*k_ins *L_p);
R_cond_pipe = R_pipe_pvc + R_pipe_ins;

% Mat per-path conduction
R_cond_top = d_mat/(k_mat*A_mat);                 % K/W
R_cond_bot = d_mat/(k_mat*A_mat);                 % K/W
R_bot_film = 1/(h_ground*A_mat);                  % K/W
R_bot      = R_cond_bot + R_bot_film;             % K/W

% ---------- Dynamic ambient ----------
rho_air = 1.2;  cp_air = 1005;  V_air_eff = 5.0;  
C_air   = rho_air*cp_air*V_air_eff;               
C_floor = 2000*800*(A_mat*0.05);                  

UA_air_to_room   = 8;   % W/K  
UA_floor_to_gnd  = 4;   % W/K  

% ---------- TIME GRID ----------
dt = 1; t = 0:dt:3730;

% ---------- STATES ----------
T_tank  = zeros(size(t));  T_tank(1)=T_ini;
T_air   = zeros(size(t));  T_air(1) = T_env0;     
T_floor = zeros(size(t));  T_floor(1)= T_env0;    
T_col_in  = zeros(size(t)); T_col_out = zeros(size(t));

% ---------- FUNCTIONS ----------
h_rad_func = @(Tsurf, eps_val, Tamb) 4*eps_val*sigmaSB*((Tsurf + Tamb)/2).^3;

R_wall_fun  = @(Tsurf, Tair) R_cond_wall + 1/((h_conv_air + h_rad_func(Tsurf,0.85,Tair))*A_wall_outer);
R_floor_fun = @(Tsurf, Tair) R_cond_floor+ 1/((h_conv_air + h_rad_func(Tsurf,0.85,Tair))*A_floor);
R_lid_fun   = @(Tsurf, Tair) R_cond_lid  + 1/((h_conv_air + h_rad_func(Tsurf,0.85,Tair))*A_lid);
R_pipe_fun  = @(Tsurf, Tair) R_cond_pipe + 1/((h_conv_air + h_rad_func(Tsurf,eps_pipe,Tair))*A_pipe_ext);

% ---------- IAM SETTINGS ----------
fb = 0.20; fd = 0.80;     % 20% beam, 80% diffuse
eta0_b = eta0; eta0_d = eta0;
b0 = 0; theta = 0;        % beam slope 0, lamp normal
Kb = 1; Kd = 1;           % IAM factors

% ---------- SOLVER ----------
for k = 1:length(t)-1
    Tair   = T_air(k);
    Tfloor = T_floor(k);

    % Mat losses
    Tmat = T_tank(k);
    h_rad_top = h_rad_func(Tmat, eps_mat, Tair);
    R_top_film = 1/((h_conv_air + h_rad_top)*A_mat);
    R_top = R_cond_top + R_top_film;

    Q_top = max((Tmat - Tair)  / R_top,  0);
    Q_bot = max((Tmat - Tfloor)/ R_bot,  0);

    % IAM-based absorbed power
    Gb = fb * G; Gd = fd * G;
    P_abs = A_mat * ( Gb * eta0_b * Kb + Gd * eta0_d * Kd );

    P_use = max(P_abs - (Q_top + Q_bot), 0);

    % Collector outlet
    dT_col = P_use / (mdot*cp);
    T_col_in(k)  = T_tank(k);
    T_col_out(k) = min(T_col_in(k) + dT_col, 372.15);

    % Tank & pipes losses
    Rw = R_wall_fun (T_tank(k), Tair);
    Rf = R_floor_fun(T_tank(k), Tair);
    Rl = R_lid_fun  (T_tank(k), Tair);
    Rp = R_pipe_fun (T_tank(k), Tair);

    Q_loss_tank = (T_tank(k) - Tair) * (1/Rw + 1/Rf + 1/Rl);
    Q_loss_pipe = (T_tank(k) - Tair) / Rp;

    % Tank ODE
    dTdt_tank = ( mdot*cp*(T_col_out(k) - T_tank(k)) - (Q_loss_tank + Q_loss_pipe) ) / (m_wt*cp);
    T_tank(k+1) = min(T_tank(k) + dTdt_tank*dt, 372.15);

    % Air node
    Q_to_air = Q_top + Q_loss_tank + Q_loss_pipe;
    dTdt_air = ( Q_to_air - UA_air_to_room*(Tair - T_env0) ) / C_air;
    T_air(k+1) = Tair + dTdt_air*dt;

    % Floor node
    dTdt_floor = ( Q_bot - UA_floor_to_gnd*(Tfloor - T_env0) ) / C_floor;
    T_floor(k+1) = Tfloor + dTdt_floor*dt;
end

% Last-point collector outlet
T_col_in(end)  = T_tank(end);
Tair   = T_air(end); Tfloor = T_floor(end); Tmat = T_tank(end);
h_rad_top = h_rad_func(Tmat, eps_mat, Tair);
R_top_film = 1/((h_conv_air + h_rad_top)*A_mat);
R_top = R_cond_top + R_top_film;
Q_top = max((Tmat - Tair)  / R_top,  0);
Q_bot = max((Tmat - Tfloor)/ R_bot,  0);

Gb = fb * G; Gd = fd * G;
P_abs = A_mat * ( Gb * eta0_b * Kb + Gd * eta0_d * Kd );

P_use = max(P_abs - (Q_top + Q_bot), 0);
T_col_out(end) = min(T_col_in(end) + P_use/(mdot*cp), 372.15);

% ---------- PLOT ----------
figure; hold on; grid on
plot(t, T_tank - 273.15, 'LineWidth',1.8, 'DisplayName','Water temp (IAM included)');
plot(t, T_air  - 273.15, '--', 'LineWidth',1.2, 'DisplayName','Air near mat');
plot(t, T_floor- 273.15, ':',  'LineWidth',1.2, 'DisplayName','Floor node');
plot(t_meas, Tin_meas,  'o-', 'DisplayName','Measured T_{in}');
plot(t_meas, Tout_meas, 's-', 'DisplayName','Measured T_{out}');
xlabel('Time [s]'); ylabel('Temperature [°C]');
xlim([0, t(end)]); legend('Location','best');
title('Model with IAM (20% beam, 80% diffuse, K_b=K_d=1)');














