function Gu
    clear;
    close all;
    clc;
    
    global Q R vr wr vm wm alpha beta
    pi = 3.1416;
    Q = diag([0.3 0.3 0.8]);
    R = diag([0.3 0.1]);
    vr = 0.4;
    wr = 0.5;
    vm = 0.5;
    wm = pi/2;
    alpha = 2;
    beta = 1;
    
    mpciterations = 50;
    N             = 10;
    T             = 0.1;
    tmeasure      = 0.0;
    xmeasure      = [ -0.2, 0, 0];
    u0            = 0.01*ones(2,N);
    tol_opt       = 1e-8;
    opt_option    = 0;
    iprint        = 5;
    type          = 'differential equation';
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;
    
    nmpc(@runningcosts, @terminalcosts, @constraints, ...
     @terminalconstraints, @linearconstraints, @system_ct, ...
     mpciterations, N, T, tmeasure, xmeasure, u0, ...
     tol_opt, opt_option, ...
     type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
     iprint, @printHeader, @printClosedloopData, @plotTrajectories);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)
    global Q R
    cost = x * Q * x' + u' * R * u;
end

function cost = terminalcosts(t, x)
    cost = 0.5 * x * x';
end

function [c,ceq] = constraints(t, x, u)
    c = [];
    ceq  = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    global vr wr vm wm
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb(1) = vr * cos(x(3)) - vm;
    ub(1) = vr * cos(x(3));
    lb(2) = wr - wm;
    ub(2) = wr + wm;
end

function [c,ceq] = terminalconstraints(t, x)
    global alpha beta vr wr vm wm
    c(1) = abs(x(2)) - abs(x(1));
    c(2) = x(2) * x(3);
    c(3) = 1/alpha * (vr * cos(x(3)) - vm) - x(1);
    c(4) = x(1) - 1/alpha * ( vr * cos(x(3)) );
    c(5) = -1/beta * (wm + wr) - x(3);
    c(6) = x(3) - 1/beta * (wm - wr);
    ceq  = [];
end

function dx = system_ct(t, x, u, T)
    global vr wr
    dx    = zeros(3,1);
    w = wr - u(2);
    dx(1) = -w * x(2) + u(1);
    dx(2) = w * x(1) + vr * sin(x(3));
    dx(3) = u(2);
end

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
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f', mpciter, ...
            u(1), u(2), x(1), x(2), x(3), t_Elapsed);
end

function plotTrajectories(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
	
	transR_inv = @(theta)[cos(theta) -sin(theta); sin(theta) cos(theta)];
	xr = 0.8 * cos(0.5*t_intermediate);
    yr = 0.8 * sin(0.5*t_intermediate);
    thetar = 0.5 * t_intermediate;
    for i = 1:1:size(t_intermediate)
        if thetar(i) >= pi
            thetar(i) = thetar(i) - pi;
        end
        theta = thetar(i) - x_intermediate(i,3);
        X(:,i) = [xr(i) yr(i)]' - transR_inv(theta) * x_intermediate(i, 1:2)';
    end
    
    figure(1);
	xlabel('t');
    hold on;
    plot(t_intermediate, x_intermediate(:,1), '-or');
    plot(t_intermediate, x_intermediate(:,2), '-og', t_intermediate, x_intermediate(:,3), '-ob');
    grid on;
    legend('xe', 'ye', '��e');
    axis([0 t_intermediate(end) -0.5 0.5]);
    
    figure(2);
    xlabel('x');
    ylabel('y');
    hold on;
    plot(xr, yr, '-.b', 'LineWidth', 2);
    plot(X(1,:), X(2,:), '-or');
    grid on;
    legend('reference', 'real');
    axis equal;
end