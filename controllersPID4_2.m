%% load following

% Declare global variable U for logging
global U U_t;

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

%dt       = 1e-4;    % ideally it should be the time elapsed between the current and previous time steps

% Initial state
x0       = [1         1         1         1         1         1         1 62803189247020.48 1018930579495656.25        900.42       898.28        888.261];
%x0       = [1          1          1          1          1          1          1 23549641141.3791 16604965156.7948        900.42        898.28        882.61]; %-0.000221924659921813];
%x0       = [0.99          0.99          0.99          0.99          0.99          0.99          0.99 62803189247020.45 1018930579495658.62        900.42        898.28        882.61]; %-0.000221924659921813];
%x0 = [0.00         -0.00         -0.00         -0.00         -0.00         -0.00         -0.00 62803189247020.48 1018930579495656.25        864.00        864.00        864.00];
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

% Certainly! This code appears to be part of a MATLAB script that implements a numerical
% integration and control algorithm for a system modeled by differential 
% algebraic equations (DAEs). Here's an explanation for each line:
% 
% 1. `x(1, :) = x0;`
%    - This line initializes the first row of the matrix `x` with the 
% initial state vector `x0`. The colon `:` operator is used to assign values 
% to all columns of the first row, indicating that `x` may be a matrix 
% where each row represents the state at a different time point.
% 
% 2. `u(1) = u0;`
%    - This initializes the control input vector `u` with the initial control value `u0`. 
% The first element of the vector `u` is set to `u0`, which is the control input 
% at the initial time step.
% 
% 3. `for i=1:length(time)-1`
%    - This line starts a `for` loop that iterates from `1` to `length(time)-1`. 
% Here, `time` is presumably a vector containing the simulation time points. 
% The loop will perform operations for each time step, except the last one.
% 
% 4. `[dx, algebraic] = reactorDAE(0, x(i,:), u(i));`
%    - This line calls a function `reactorDAE` with inputs `0`, the current state `x(i,:)`, 
% and the current control input `u(i)`. The function returns `dx`, the derivative of the state, 
% and `algebraic`, some algebraic variables which are not used further in this snippet. 
% The `0` might represent the current time (if the function needs it), or 
% it could be a placeholder for another parameter.
% 
% 5. `x(i+1,:) = x(i,:) + dx'*dt;`
%    - Here, the next state `x(i+1,:)` is computed using Euler's method for numerical integration.
% `dx'` is the transpose of `dx`, ensuring correct dimensionality (if needed), 
% and `dt` is the time step size. The current state `x(i,:)` is updated by adding the 
% product of the derivative `dx` and the time step `dt`.
% 
% 6. `u(i+1) = pidController(time(i), x(i+1,1), ref(i+1), Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 0, 0);`
%    - This line computes the next control input `u(i+1)` by calling 
% a PID controller function `pidController`. The inputs to this function are the current time 
% `time(i)`, the next state value `x(i+1,1)` (presumably the first component of the state vector),
% the reference signal value `ref(i+1)`, PID gains `Kp`, `Ki`, `Kd`, anti-windup parameter `Kaw`, 
% a time constant `T_C`, limits `max`, `min`, the maximum rate of change `max_rate`,
% and two zeros which might be placeholders for additional parameters.
% 
% 7. `t = time;`
%    - This line simply assigns the vector `time` to the variable `t`, 
% possibly for use outside the loop or for clarity and simplicity in later code.
% 
% 8. `du = [0 diff(u)/dt];`
%    - Finally, this line computes the derivative of the control input `u` with respect to time.
% `diff(u)` computes the differences between consecutive elements of `u`,
% which approximates the derivative, and dividing by `dt` scales it by the time step. 
% The vector is prepended with a `0` to maintain the same length as `u`.
% 
% This code snippet is typical for simulations of dynamical systems where 
% states are updated based on differential equations and controlled 
% by a feedback mechanism such as a PID controller.

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
    % nv    = 2.2e+3;           % average velocisty of neutrons (thermal)
    % nu    = 2.43;             % the total number of neutrons liberated per rx
    % Sigf  = 1/(gen*nv*nu);    % what is this 
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

%     degrees = 0:10:350; % Control drum positions from 0 to 350 degrees, every 10 degrees
% 
%     reactivityPerTenDegrees = [
%     -0.033085599
%     -0.032318805
%     -0.031579751
%     -0.029910666
%     -0.027954186
%     -0.0253718
%     -0.022016275
%     -0.018495885
%     -0.014448794
%     -0.010209174
%     -0.005749872
%     -0.001313724
%     0.002652943
%     0.006467894
%     0.009567577
%     0.01191629
%     0.013242288
%     0.01364134
%     0.013281235
%     0.011535382
%     0.009547958
%     0.006369174
%     0.002722567
%     -0.001473167
%     -0.005864188
%     -0.010443949
%     -0.014540393
%     -0.018772929
%     -0.022302552
%     -0.025346567
%     -0.027783031
%     -0.029751583
%     -0.031356326
%     -0.032480814
%     -0.033317248
%     -0.033157111
% ];

    %Gr4     = 1450E-5/180*2; % half
    % Gr2   = -660E-5/180*2; % 1/4
    %Gr1     = 250E-5/180*2; % one
    Gr4=26.11e-5;
    % angles=10:10:360;
    % current_angle=180;
    % Gr=interp1(angles,G_V,current_angle,'linear');
    %Gr=1450E-5/180*2;
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
    %Rho_d1  = x(13);
    %Rho_d2  = x(14);
    
    % u is the control input = position of the control drum in degrees - interpolate to find the reactivity
    %Rho_d1 = interp1(degrees, reactivityPerTenDegrees , u,'linear');
    Rho_d0=-0.033085599;
    Rho_d1 =Rho_d0+u*Gr4;
    
    G       = 3.2e-11;
    V       = 400*200; 
    Pi      = P_0/(G*Sum_f*V);
    
    %% The extra parameter u in reactorDAE(t, x, u) comes from the pidController()
    %u1  = u;    % actuation channel 
    %u4  = u;    % actuation thru channel x14
    
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

    % Store control input globally (for plotting)
    %U = [U, u];
    %U_t = [U_t, t];

end
