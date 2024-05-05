%% load following


dt=0.1;
T = 2000;
n0=1;
time=0:dt:T;
nt=length(time);
ref=zeros(nt,1);

%time_point=[0 20 30 50 60 80 90 110 130 150];%work with only this time
%time_point=[0 10 40 70 80 100 110 130 140 150];
%time_point=[0 30 40 50 60 70 90 100 140 150];
time_point=[0 20 30 50 60 80 90 110 130 200]*10;
%pow_point=[0.8 0.8 0.4 0.4 0.8 0.8 0.4 0.4 0.8 0.8];power Transient
pow_point=[1 1 0.4 0.4 1 1 0.4 0.4 1 1];
%pow_point=[0.3 1 0.4 0.4 0.6 0.6 0.8 0.8 1 1];
%pow_point=[1 1 1 1 1 1 1 1 1 1];
%pow_point=[0.9 0.9 0.3 0.3 0.9 0.9 0.3 0.3 0.9 0.9];
%pow_point=[0.7 0.7 0.4 0.4 0.4 0.4 0.6 0.6 0.7 0.7];
%pow_point=[0.5 0.5 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];

ref_old=pow_point(1);
for it=1:nt %to vectorize nt or make each point of nt
    if(it>1)
        time(it)=time(it-1)+dt;%to be able to make equal time step of the time
        ref_old=ref(it-1); %for each time slot 
    end
    ref(it)=ref_old;
    i1=1;
    i2=1;
    for ii=1:length(time_point)-1
        if(time_point(ii)<=time(it) && time(it)<=time_point(ii+1))
            i1=ii;
            i2=ii+1;
            frac1=(time_point(ii+1)-time(it))/(time_point(ii+1)-time_point(ii));
            frac2=1.0-frac1;
            ref(it)=frac1*pow_point(i1)+frac2*pow_point(i2);
            break
        end
    end

end
ref(:)=ref(:);
yref=(ref.*n0);

setpoint = yref;

% PID tuning - start small
Kp          = 1;          % Proportional gain
Ki          = 1;           % Integral gain
Kaw         = 0.3;          % Anti-windup gain
Kd          = 0.001;           % Derivative gain
T_C         = 0.2;
max         = 180;
min         = 0;
max_rate    = 0.5;
%max_rate    = 1;


% Initial state
x0       = [1         1         1         1         1         1         1 62803189247020.48 1018930579495656.25        900.42       898.28        888.261];
     
% Initialise PID controller internal states
%u0 = 83.84;
u0=77.5556;
pidController(0, 0, 0, Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 1, u0);

% Call ode45
%[t, x] = ode45(@(t, x) reactorDAE(t, x, pidController(t, x(1), setpoint(round(t/dt)+1), Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 0, 0)), [0 time(end)], x0);
%[t, x] = ode15s(@(t, x) reactorDAE(t, x, pidController(t, x(1), setpoint(round(t/dt)+1), Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 0)), [0 time(end)], x0);
%[t, x] = ode23s(@(t, x) reactorDAE(t, x, pidController(t, x(1), setpoint(round(t/dt)+1), Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 0)), [0 time(end)], x0);


%% Simulation
x(1, :) = x0;
u(1) = u0;
for i=1:length(time)-1
    [dx, algebraic] = reactorDAE(0, x(i,:), u(i));
    x(i+1,:) = x(i,:) + dx'*dt;

    u(i+1) = pidController(time(i), x(i+1,1), ref(i+1), Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 0, 0);
end
t = time;
du = [0 diff(u)/dt];



%% Interpolate ref in case t and time are different
ref_interpolated=interp1(time,ref,t,'linear');

%% Plot result
figure(1)

% Power subplot
subplot(3, 1, 1);
plot(t, x(:,1)*100, 'LineWidth', 2);
hold on;
plot(t,ref_interpolated*100, '--', 'LineWidth', 2);
hold off;
grid on;
xlabel('Time (s)');
ylabel('Power (%)');
title('PID microcontroller system core power simulation');

ylim([0,200]);
legend({'Actual power', 'Desired power'});

% Rotation subplot
subplot(3, 1, 2);
plot(time, u, 'LineWidth', 2), grid on, xlabel('Time (s)'), ylabel('Rotation (deg)'), title('Rotation');
ylim([70,80])

% Rotation rate of change subplot
subplot(3, 1, 3);
plot(time, du, 'LineWidth', 2), grid on, xlabel('Time (s)'), ylabel('Rate of change (deg/s)'), title('Rate of change');
ylim([-0.5,0.5])
saveas(figure(1),'power_simulation.png')
%% Reactor dynamics
function [dx, algebraic] = reactorDAE(t, x, u)
    %% Parameters
    Sig_x   = 2.65e-22;
    yi      = 0.061;
    yx      = 0.002;
    lamda_x = 2.09e-5;
    lamda_I = 2.87e-5;
    Sum_f   = 0.3358;
    
    l       = 1.68e-3;
    beta    = 0.0048;
    beta_1  = 1.42481E-04; 
    beta_2  = 9.24281E-04;
    beta_3  = 7.79956E-04;
    beta_4  = 2.06583E-03;
    beta_5  = 6.71175E-04;
    beta_6  = 2.17806E-04;
    Lamda_1 = 1.272E-02;
    Lamda_2 = 3.174E-02;
    Lamda_3 = 1.160E-01;
    Lamda_4 = 3.110E-01;
    Lamda_5 = 1.400E+00;
    Lamda_6 = 3.870E+00;


    Reactivity_per_degree=26.11e-5;
   
    cp_f    = 977;
    cp_m    = 1697;
    cp_c    = 5188.6;
    M_f     = 2002;
    M_m     = 11573;
    M_c     = 500;
    mu_f    = M_f*cp_f;
    mu_m    = M_m*cp_m;
    mu_c    = M_c*cp_c;
    f_f     = 0.96;
    P_0     = 22e6; 
    Tf0     = 1105;
    Tm0     = 1087;
    T_in    = 864;
    T_out   = 1106;
    Tc0     = (T_in+T_out)/2;
    K_fm    = f_f*P_0/(Tf0-Tm0);
    K_mc    = P_0/(Tm0-Tc0);
    M_dot   = 1.75E+01;
    alpha_f = -2.875e-5;
    alpha_m = -3.696e-5;
    alpha_c = 0.0;
    X0      = 2.35496411413791e10;

    %% Declaration of state variables, x(i), where i = 1 to 14
    n_r     = x(1); 
    Cr1     = x(2); 
    Cr2     = x(3); 
    Cr3     = x(4); 
    Cr4     = x(5);
    Cr5     = x(6);
    Cr6     = x(7);
    X       = x(8); 
    I       = x(9); 
    Tf      = x(10); 
    Tm      = x(11); 
    Tc      = x(12); 
    
    Rho_d0=-0.033085599;
    Rho_d1 =Rho_d0+u*Reactivity_per_degree;
    
    G       = 3.2e-11;
    V       = 400*200; 
    Pi      = P_0/(G*Sum_f*V);
    
    %% ODEs
    dx      = zeros(12,1);
    rho     = Rho_d1 + alpha_f*(Tf - Tf0) + alpha_c*(Tc - Tc0) + alpha_m*(Tm - Tm0) - Sig_x*(X - X0)/Sum_f;
    
    %% kinetics equations with six-delayed neutron groups
    dx(1)   = (rho-beta)/l*n_r+beta_1/l*Cr1+beta_2/l*Cr2+beta_3/l*Cr3+beta_4/l*Cr4+beta_5/l*Cr5+beta_6/l*Cr6;
    dx(2)   = Lamda_1*n_r-Lamda_1*Cr1;
    dx(3)   = Lamda_2*n_r-Lamda_2*Cr2;
    dx(4)   = Lamda_3*n_r-Lamda_3*Cr3;
    dx(5)   = Lamda_4*n_r-Lamda_4*Cr4;
    dx(6)   = Lamda_5*n_r-Lamda_5*Cr5;
    dx(7)   = Lamda_6*n_r-Lamda_6*Cr6;
    
    %% Xenon and Iodine dynamics
    dx(8)   = yx*Sum_f*Pi+lamda_I*I-Sig_x*X*Pi-lamda_x*X;
    dx(9)   = yi*Sum_f*Pi-lamda_I*I;
    
    %% thermal–hydraulics model of the reactor core
    %dt       = 1e-4; 
    dx(10)  = f_f*P_0/mu_f*n_r-K_fm/mu_f*(Tf-Tc);
    dx(11)  = (1-f_f)*P_0/mu_m*n_r+(K_fm*(Tf-Tm)-K_mc*(Tm-Tc))/mu_m;
    dx(12)  = K_mc*(Tm-Tc)/mu_c-2*M_dot*cp_c*(Tc-T_in)/mu_c;
    %dx(13)  = Gr(round(t/dt)+1)*u1;
    %dx(13)  = Gr*u1;
    
    %Algebraic equations
    algebraic = [];
end

%% PID controller
function u = pidController(t, measurement, setpoint, Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, init, u0)

    persistent integral err_prev deriv_prev t_prev command_sat_prev command_prev;
    %global U U_t;

    % Initialise static variables    
    if (isempty(integral) || init == 1)
        integral         = u0;
        err_prev         = 0;
        t_prev           = 0;
        deriv_prev       = 0;
        command_prev     = u0;
        command_sat_prev = u0;
        %U = [];
        %U_t = [];
    end
    
    % Calculate error
    err = setpoint - measurement;

    % Calculate period
    T = t - t_prev;

    % Store previous time
    t_prev = t;
    
    % Update integral term with anti-windup
    integral = integral + Ki*err*T + Kaw*(command_sat_prev - command_prev)*T;
    
    % Calculate filtered derivative
    deriv_filt = (err - err_prev + T_C*deriv_prev) / (T + T_C);

    % Store previous error and previous derivative
    err_prev = err;
    deriv_prev = deriv_filt;
    
    % Calculate command using PID equation
    command = Kp*err + integral + Kd*deriv_filt;
    
    % Store previous command
    command_prev = command;
    
    % Saturate command
    if command > max
        command_sat = max;
    elseif command < min
        command_sat = min;
    else
        command_sat = command;
    end
    
    % Apply rate limiter
    if command_sat > command_sat_prev + max_rate*T
        command_sat = command_sat_prev + max_rate*T;
    elseif command_sat < command_sat_prev - max_rate*T
        command_sat = command_sat_prev - max_rate*T;
    end
    
    % Store previous saturated command
    command_sat_prev = command_sat;
    
    % Output the saturated command
    u = command_sat;

end
