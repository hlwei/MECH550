%%for case1
% tracking circle shape

function UUVTracking
    clear;
    close all;
    clc;
    
    global Q R P
    Q = diag([10^5 10^5 10^3 10^2 10^2 10^2]);
    R = diag([10^(-4) 10^-4 10^-4]);
    P = diag([10^3 10^3 10^2 10 10 10]);
    
    %X = [];
    mpciterations = 300;
    N             = 5;
    T             = 0.1;
    tmeasure      = 0.0;
    xmeasure      = [-1 0 0 0 0 0];
    u0            = 0.0*ones(3,N);
    tol_opt       = 1e-8;
    opt_option    = 0;
    iprint        = 5;
    type          = 'differential equation';
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;
    
    [t,x,u]=nmpc(@runningcosts, @terminalcosts, @constraints, ...
     @terminalconstraints, @linearconstraints, @system_ct, ...
     mpciterations, N, T, tmeasure, xmeasure, u0, ...
     tol_opt, opt_option, ...
     type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
     iprint, @printHeader, @printClosedloopData, @plotTrajectories);
%      X = [X,x]
%      xd = 0.5*t;
%      yd = sin(0.5*t);
%      figure(1);
%      hold on
%      plot(x(1),x(2),'b');
%      plot(xd,yd,'k'); 
%      grid on;
%      
%      figure(2);
%      hold on
%      plot(xd,yd,'Linewidth',2,'-k')
     %input control
     figure(3);
     hold on
     plot(t,u(1,:));
     plot(t,u(2,:));
     plot(t,u(3,:));
     legend('u1','u2','u3');
     xlabel('Time[s]');
     ylabel('Control Input[N]');
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)
    global Q R 
    xd = sin(0.5*t);
    yd = sin(0.25*t);
    psid = atan(0.5*t);
    ud = sqrt((0.5*cos(0.5*t))^2+(0.25*cos(0.25*t))^2);
    vd = 0;
    rd = (0.5*0.25^2*cos(0.5*t)*sin(0.25*t)-0.25*0.5^2*sin(0.5*t)*cos(0.25*t))/ud^2;
    xe = x-[xd yd psid ud vd rd];
    cost = xe * Q * xe' + u' * R * u;
end

function cost = terminalcosts(t, x)
    global P 
    xd = sin(0.5*t);
    yd = sin(0.25*t);
    psid = atan(0.5*t);
    ud = sqrt((0.5*cos(0.5*t))^2+(0.25*cos(0.25*t))^2);
    vd = 0;
    rd = (0.5*0.25^2*cos(0.5*t)*sin(0.25*t)-0.25*0.5^2*sin(0.5*t)*cos(0.25*t))/ud^2;
    xe = x-[xd yd psid ud vd rd];
    cost = xe * P * xe';
end

function [c,ceq] = constraints(t, x, u)
    c = [];
    ceq  = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, ~)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb = [-1200 -1000 -500];
    ub = [1200 1000 500]; 
end

function [c,ceq] = terminalconstraints(t, x)
    c = [];
    ceq  = [];
end

%% model from Dongbin Gu 
% function dx = system_ct(t, x, u, T)
%     global vr wr
%     dx    = zeros(3,1);
%     w = wr - u(2);
%     dx(1) = -w * x(2) + u(1);
%     dx(2) = w * x(1) + vr * sin(x(3));
%     dx(3) = u(2);
% end

%% UUV model
function dx = system_ct(t,x,u,T)

% System coefficients
%==========================================================================
err_model=0;

m=116;
Iz=13.1;
X_udot=-167.6*(1+err_model);
Y_vdot=-477.2*(1+err_model);
N_rdot=-15.9*(1+err_model);
Xu=26.9*(1+err_model);
Yv=35.8*(1+err_model);
Nr=3.5*(1+err_model);
Du=241.3*(1+err_model);
Dv=503.8*(1+err_model);
Dr=76.9*(1+err_model);

Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;
%==========================================================================

dx  = zeros(6,1);
dx(1) = x(4)*cos(x(3)) - x(5)*sin(x(3));
dx(2) = x(4)*sin(x(3)) + x(5)*cos(x(3));
dx(3) = x(6);
dx(4) = (My/Mx)*x(5)*x(6)-(Xu/Mx)*x(4)-(Du/Mx)*x(4)*abs(x(4))+u(1)/Mx;
dx(5) = -(Mx/My)*x(4)*x(6)-(Yv/My)*x(5)-(Dv/My)*x(5)*abs(x(5))+u(2)/My;
dx(6) = ((Mx-My)/Mpsi)*x(4)*x(5)-(Nr/Mpsi)*x(6)-(Dr/Mpsi)*x(6)*abs(x(6))+u(3)/Mpsi;

% x_dot=u*cos(psi)-v*sin(psi);
% y_dot=u*sin(psi)+v*cos(psi);
% psi_dot=r;
% u_dot=(My/Mx)*v*r-(Xu/Mx)*u-(Du/Mx)*u*abs(u)+Fu/Mx;
% v_dot=-(Mx/My)*u*r-(Yv/My)*v-(Dv/My)*v*abs(v)+Fv/My;
% r_dot=((Mx-My)/Mpsi)*u*v-(Nr/Mpsi)*r-(Dr/Mpsi)*r*abs(r)+Fr/Mpsi;

% reference model

% global ud vd rd
% [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
%                                           x0, u, atol_ode, rtol_ode, type);
% 	dx    = zeros(6,1);
% 	xd = 1 * sin(0.5*t_intermediate);
%     yd = 1 * sin(0.25*t_intermediate);
%     psid = atan2(0.5 * t_intermediate);
%     ud = 1
%     vd = 0
%     rd = 0.8
% kinematics 
% psi = x(6)+psid
% Fud = Mx*dotud-My*vd*rd+Du*ud
% Frd = Mpsi*dotrd+(My-Mx)*ud*vd+Dr*rd
% 
% dx(1) = x(4)*cos(psi)-sin(psi)*x(5)+(cos(psi)-cos(psid))*ud+(-sin(psi)+sin(psid))*vd; 
% dx(2) = sin(psi)*x(4)+cos(psi)*x(5)+(sin(psi)-sin(psid))*ud+(cos(psi)-cos(psid))*vd;
% dx(3) = x(6);
% dx(4) = My/Mx*(x(5)*x(6)+x(5)*rd+vd*x(6))-(Du/Mx)*x(4)+(u(1)-Fud)/Mx;
% dx(5) = -Mx/My*(x(4)*x(6)+x(4)*rd+ud*x(6));
% dx(6) = (Mx-My)/Mpsi*(x(4)*x(5)+x(4)*vd+ud*x(5))-(Dr/Mpsi)*x(6)+(Fr-Frd)/Mpsi;
%X_dot=[x_dot;y_dot;psi_dot;u_dot;v_dot;r_dot];
%dx = X_dot;
%Xplus=X+X_dot*T;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of output format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printHeader()
    fprintf('   k  |      u1(k)        u2(k)        x(1)        x(2)        x(3)        Time\n');
    fprintf('-------------------------------------------------------------------------------\n');
end

function printClosedloopData(mpciter, u, x, t_Elapsed)
    global t_total 
    t_total = t_total + t_Elapsed;
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f', mpciter, ...
            u(1), u(2),u(3), x(1), x(2), x(3), t_Elapsed);
end

function plotTrajectories(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
	
	%transR_inv = @(theta)[cos(theta) -sin(theta); sin(theta) cos(theta)];
    %the reference trajectory---eight shape
%     xr = 1 * sin(0.5*t_intermediate);
%     yr = 1 * sin(0.25*t_intermediate);
%     thetar = atan(0.5 * t_intermediate);
      xr = sin(0.5*t_intermediate);
      yr = sin(0.25*t_intermediate);
      thetar = atan(0.5*t_intermediate);
     for i = 1:1:size(t_intermediate)
        if thetar(i) >= pi
            thetar(i) = thetar(i) - pi;
        end
          thetae(i) = thetar(i) - x_intermediate(i,3);
          Xe(i,:) = [xr(i) yr(i)] -  x_intermediate(i, 1:2);
     end
    %size(Xe)
    %size(thetae)
    figure(1);
    xlabel('Time[s]');
    ylabel('Error');
    hold on;
    %plot(x(1),x(2),'k*');
    plot(t_intermediate, Xe(:,1), '-r','LineWidth', 1.5);
    plot(t_intermediate, Xe(:,2), '-g','LineWidth', 1.5);
    plot(t_intermediate, thetae, '-b','LineWidth', 1.5);
    grid on;
    legend('x_e[m]', 'y_e[m]', '\theta_e[rad]');
    axis([0 t_intermediate(end) -2 2]);
    %size(xr)
    %size(x_intermediate)
    
    figure(2);
    xlabel('x[m]');
    ylabel('y[m]');
    hold on;
    plot(xr, yr, '-k', 'LineWidth', 1.5);
    plot(x_intermediate(:,1), x_intermediate(:,2), '-r', 'LineWidth', 1.5);
    grid on;
    legend('reference trajectory', 'real trajectory');
    axis equal;

    
end
