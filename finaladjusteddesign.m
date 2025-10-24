%% ===================== PART A — Read group05.txt (Tin only, seconds) ====
fname = 'group05.txt';
fid = fopen(fname, 'r'); assert(fid > 0, 'Could not open %s', fname);
C = textscan(fid, '%s %f %f %f %f', ...
    'HeaderLines', 1, 'Delimiter', {' ', '\t', ';'}, ...
    'MultipleDelimsAsOne', true, 'CollectOutput', true);
fclose(fid);

timeStr   = C{1};
nums      = C{2};
Tin_txt   = nums(:,1);    % use Tin only

% Convert HH:mm:ss to elapsed seconds from first sample
t_txt_dt  = datetime(timeStr, 'InputFormat', 'HH:mm:ss');
t0_dt     = t_txt_dt(1);
t_txt_sec = seconds(t_txt_dt - t0_dt);

%% ===================== PART B — CPC + Greenhouse model ==================
rng(12,'twister');

% --- OPTIONAL: shift model time to align with data (e.g., pump-on delay)
offset_seconds = 0;   % set +/− seconds as needed

% --------- CPC geometry / properties (geometry known/unchanged) ----------
D_in = 0.600; D_out = 0.400; L = 0.20;
Rin = D_in/2; Rout = D_out/2;
A_in = pi*Rin^2; A_out = pi*Rout^2;

% Optics (separated optical vs thermal)
rho         = 0.85;         % CPC wall reflectance
I_sun_total = 1000;         % W/m^2
tau_glass   = 0.92;         % glazing transmittance
alpha_absorb = 0.97;        % tube solar absorptance
eps_emit     = 0.80;        % tube thermal emissivity
solar_gain_factor = 1.10;   % mild global gain correction (+10%)

% Parabola for CPC profile
delta = Rin - Rout;  a_par = -delta/(L^2);  b_par = 2*delta/L;
r_fun = @(z) a_par.*(z.^2) + b_par.*z + Rout;
drdz  = @(z) 2*a_par.*z + b_par;

% Monte Carlo settings (ray trace to estimate CPC throughput)
N_total = 20000; frac_parallel = 0.20;
N_par = round(frac_parallel*N_total); N_diff = N_total - N_par;
max_reflections = 50; energy_tol = 1e-6;
Rsrc = 1.5*Rin; theta_max = pi/2;
sample_disk = @(N,R) sqrt(rand(N,1)).*R.*exp(1i*2*pi*rand(N,1));

sum_exit = 0;
% Parallel rays
pts = sample_disk(N_par, Rsrc);
for i=1:N_par
    pos=[real(pts(i)); imag(pts(i)); L]; dir=[0;0;-1];
    [w_out, ~, ~] = trace_ray(pos, dir, 1, L, Rin, Rout, r_fun, drdz, rho, energy_tol, max_reflections);
    sum_exit = sum_exit + w_out;
end
% Diffuse rays
pts = sample_disk(N_diff, Rsrc);
u=rand(N_diff,1);
theta_d = acos(1 - u*(1 - cos(theta_max)));
phi_d = 2*pi*rand(N_diff,1);
ux = sin(theta_d).*cos(phi_d);
uy = sin(theta_d).*sin(phi_d);
uz = -cos(theta_d);
for i=1:N_diff
    pos=[real(pts(i)); imag(pts(i)); L]; dir=[ux(i); uy(i); uz(i)];
    [w_out, ~, ~] = trace_ray(pos, dir, 1, L, Rin, Rout, r_fun, drdz, rho, energy_tol, max_reflections);
    sum_exit = sum_exit + w_out;
end
efficiency = sum_exit / N_total;
P_in  = I_sun_total * A_in;
P_out = efficiency * P_in;        % total power through exit
E_exit = P_out / A_out;           % exit-plane irradiance (W/m^2)

fprintf('CPC Efficiency = %.4f\n', efficiency);
fprintf('P_out (total)  = %.2f W (incident %.2f W)\n', P_out, P_in);

% --------------------- Greenhouse / tank model --------------------------
% Flow (reduced to increase ΔT and linearity)
Vdot = 0.75 * (1.48/60*1e-3);   % m^3/s (25% lower than spec)

% Ambient and initial water temperature
T_env0 = 22.7 + 273.15;          % K
T_ini  = 22.7 + 273.15;        % K

% Copper tube geometry (fixed)
L_tube = 5.0; D_tube_o = 0.012; D_tube_i = 0.010;
ro = D_tube_o/2; ri = D_tube_i/2;
A_Coppertube = pi * D_tube_o * L_tube;
A_tube_proj  = D_tube_o * L_tube;

% Tank / materials
V_wt=2e-3; H_t=0.20; D_t=0.11; d_wall=0.003; d_lid=0.005; d_1=0.05;
k_1=0.04; k_pvc=0.19; rho_w=0.997e3; cp=4184;

% Loss parameters (further reduced to avoid end flattening)
h_conv_air   = 2.0;   % W/m^2K
h_ground     = 70;    % W/m^2K
contact_fraction = 0.40;
UA_air_to_room  = 2.0;  % W/K
UA_floor_to_gnd = 2.0;  % W/K

% Pipe insulation (kept good)
L_p=3.0; Di_p=0.010; t_p=0.0015;
t_ins_p = 0.020;    % m
k_pipe  = 0.19;     % W/mK (PVC)
k_ins   = 0.035;    % W/mK
eps_pipe=0.85;

% Constants
sigmaSB = 5.670374419e-8;
rho_air = 1.2; cp_air = 1005;

% Derived
mdot = rho_w*Vdot;  m_wt = rho_w*V_wt;
r1=D_t/2; r2=r1+d_wall; r3=r2+d_1;
A_floor = pi*r1^2; A_wall_outer = 2*pi*r3*H_t; A_lid = A_floor;

R_wall_pvc = log(r2/r1)/(2*pi*k_pvc*H_t);
R_wall_ins = log(r3/r2)/(2*pi*k_1*H_t);
R_cond_wall = R_wall_pvc + R_wall_ins;
R_cond_floor = d_wall/(k_pvc*A_floor) + d_1/(k_1*A_floor);
R_cond_lid   = d_lid /(k_pvc*A_lid);

% Pipe conduction
ri_p=Di_p/2; ro_p=ri_p+t_p; ro_i=ro_p+t_ins_p; A_pipe_ext=2*pi*ro_i*L_p;
R_pipe_pvc = log(ro_p/ri_p)/(2*pi*k_pipe*L_p);
R_pipe_ins = log(ro_i/ro_p)/(2*pi*k_ins*L_p);
R_cond_pipe = R_pipe_pvc + R_pipe_ins;

% Copper tube conduction (cylindrical shell)
k_copperTube = 400;
R_cond_tube_wall = log(ro/ri) / (2*pi*k_copperTube*L_tube);
R_cond_top = R_cond_tube_wall;   % top conduction through tube wall
R_cond_bot = R_cond_tube_wall;

% Bottom film/contact
A_contact = (L_tube*pi*D_tube_o)/2 * contact_fraction;
R_bot_film = 1 / (h_ground * A_contact);
R_bot = R_cond_bot + R_bot_film;

% Air/floor capacities
V_air_eff = 0.4 * 0.5 * 1.7;     % m^3 (smaller inertia)
C_air   = rho_air*cp_air*V_air_eff;
C_floor = 2000*800*(A_Coppertube*0.05);  % crude slab capacity

% Time grid (seconds)
dt = 1;
t_sim = 0:dt:round(max(t_txt_sec));   % simulate through data span
t_sim = t_sim + offset_seconds;

% Radiation helper
h_rad_func = @(Tsurf, eps_val, Tamb) 4*eps_val*sigmaSB*((Tsurf + Tamb)/2).^3;

R_wall_fun  = @(Tsurf, Tair) R_cond_wall + 1/((h_conv_air + h_rad_func(Tsurf,0.85,Tair))*A_wall_outer);
R_floor_fun = @(Tsurf, Tair) R_cond_floor+ 1/((h_conv_air + h_rad_func(Tsurf,0.85,Tair))*A_floor);
R_lid_fun   = @(Tsurf, Tair) R_cond_lid  + 1/((h_conv_air + h_rad_func(Tsurf,0.85,Tair))*A_lid);
R_pipe_fun  = @(Tsurf, Tair) R_cond_pipe + 1/((h_conv_air + h_rad_func(Tsurf,eps_pipe,Tair))*A_pipe_ext);

% States
T_tank  = zeros(size(t_sim));  T_tank(1)=T_ini;
T_air   = zeros(size(t_sim));  T_air(1) = T_env0;
T_floor = zeros(size(t_sim));  T_floor(1)=T_env0;

% Glazing opening area (geometry fixed)
A_glass_open = 0.4^2*pi/4;
frac_air = 0.10;   % fraction of glass gain straight to air

for k = 1:length(t_sim)-1
    Tair   = T_air(k);
    Tfloor = T_floor(k);

    % Top convection + radiation for tube surface (use eps_emit)
    h_rad_top = h_rad_func(T_tank(k), eps_emit, Tair);
    R_top_film = 1/((h_conv_air + h_rad_top)*A_Coppertube);
    R_top = R_cond_top + R_top_film;

    % Heat flows from tank/tube
    Q_top = (T_tank(k) - Tair)   / R_top;
    Q_bot = (T_tank(k) - Tfloor) / R_bot;

    % Solar absorbed by tube (with gain factor)
    P_abs_from_CPC   = solar_gain_factor * (E_exit * A_tube_proj * alpha_absorb);
    P_total_on_glass = solar_gain_factor * (tau_glass * I_sun_total * A_glass_open);
    P_additional_to_tube = (1-frac_air) * P_total_on_glass;
    P_abs = P_abs_from_CPC + P_additional_to_tube;

    % Net to water after losses to air/floor
    P_use = P_abs - (Q_top + Q_bot);

    % Collector outlet (water temperature rise)
    dT_col = P_use / (mdot*cp);
    T_col_out = min(T_tank(k) + dT_col, 372.15);

    % Tank & pipes losses
    Rw = R_wall_fun(T_tank(k), Tair);
    Rf = R_floor_fun(T_tank(k), Tair);
    Rl = R_lid_fun(T_tank(k), Tair);
    Rp = R_pipe_fun(T_tank(k), Tair);
    Q_loss_tank = (T_tank(k) - Tair) * (1/Rw + 1/Rf + 1/Rl);
    Q_loss_pipe = (T_tank(k) - Tair) / Rp;

    % Tank ODE
    dTdt_tank = ( mdot*cp*(T_col_out - T_tank(k)) - (Q_loss_tank + Q_loss_pipe) ) / (m_wt*cp);
    T_tank(k+1) = min(T_tank(k) + dTdt_tank*dt, 372.15);

    % Greenhouse air node
    P_air_solar = frac_air * P_total_on_glass;   % W directly into air
    Q_to_air = Q_top + Q_loss_tank + Q_loss_pipe;
    dTdt_air = (Q_to_air + P_air_solar - UA_air_to_room*(Tair - T_env0) ) / C_air;
    T_air(k+1) = Tair + dTdt_air*dt;

    % Floor node
    dTdt_floor = ( Q_bot - UA_floor_to_gnd*(Tfloor - T_env0) ) / C_floor;
    T_floor(k+1) = Tfloor + dTdt_floor*dt;
end

fprintf('Final water temperature: %.2f °C\n', T_tank(end) - 273.15);

%% ===================== PART C — Single plot (seconds) ===================
figure('Color','w'); hold on; grid on; box on;
plot(t_txt_sec, Tin_txt, 'LineWidth', 1.8, 'DisplayName', 'Practical results');
plot(t_sim, T_tank - 273.15, '--', 'LineWidth', 2.2, 'DisplayName', 'Model temp');
xlabel('Time (s)'); ylabel('Temperature (°C)');
title('Test vs Model Tank Temperature');
legend('Location','best'); set(gca,'FontSize',11);

%% ===================== LOCAL FUNCTION ===============================
function [w_exit, path, success] = trace_ray(pos,d,w,L,Rin,Rout,r_fun,drdz,rho,energy_tol,max_reflections)
    path = pos'; w_exit = 0; reflections = 0; success = false;
    while true
        if abs(d(3)) < 1e-14, t_exit = inf; else, t_exit = (0 - pos(3)) / d(3); end
        if t_exit > 1e-12
            x_exit = pos(1) + d(1)*t_exit; y_exit = pos(2) + d(2)*t_exit;
            if (x_exit^2 + y_exit^2) <= Rout^2
                w_exit = w; path(end+1,:) = [x_exit y_exit 0]; success = true; return
            end
        end
        fun = @(t) sqrt((pos(1)+d(1)*t).^2 + (pos(2)+d(2)*t).^2) - r_fun(pos(3)+d(3)*t);
        t_min=1e-8; t_max=5*L; ts=linspace(t_min,t_max,200); vals=arrayfun(fun,ts); br=[];
        for k = 1:numel(vals)-1
            if (vals(k)*vals(k+1) <= 0) && isfinite(vals(k)) && isfinite(vals(k+1))
                br = [ts(k), ts(k+1)]; break
            end
        end
        if isempty(br), return; end
        try t_wall = fzero(fun, br); catch, return; end
        if t_wall <= 1e-10, return; end
        pos = pos + d*t_wall; path(end+1,:) = pos';
        rho_r = norm(pos(1:2)); if rho_r < 1e-12, return; end
        dzdr = drdz(pos(3)); n = [pos(1)/rho_r; pos(2)/rho_r; -dzdr]; n = n/norm(n);
        d = d - 2*(d'*n)*n; w = w*rho; reflections = reflections+1;
        if w < energy_tol || reflections > max_reflections, return; end
        if d(3) > 0
            t_ent = (L - pos(3)) / d(3);
            if t_ent > 1e-12
                x_ent = pos(1) + d(1)*t_ent; y_ent = pos(2) + d(2)*t_ent;
                if (x_ent^2 + y_ent^2) <= Rin^2
                    pos = [x_ent; y_ent; L]; path(end+1,:) = pos'; return
                end
            end
        end
    end
end
