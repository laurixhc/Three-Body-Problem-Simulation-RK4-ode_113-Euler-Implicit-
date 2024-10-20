global G m1 m2 m3
 G = 1;
% Define masses
m1 = 1;
m2 = 1;
m3 = 1;

% Define initial positions and velocities
p1_start = [0.97000436, -0.24308753];
v1_start = [0.93240737/2, 0.86473146/2];

p2_start = [-0.97000436, 0.24308753];
v2_start = [0.93240737/2, 0.86473146/2];

p3_start = [0,0];
v3_start = [-0.93240737, -0.86473146];

% Define the time range
t_initial = 0.0;
t_final = 4.0;
t_range = [t_initial, t_final];

% Concatenate initial positions and velocities into a single vector
x_initial = [p1_start, p2_start, p3_start, v1_start, v2_start, v3_start];

%-------------------------------------------------------------------------------------------------------------------
% Matlab's ode_113

% Set error tolerances for ode113
options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0E-10);

% Call ode113 to solve the system of ODEs
[T, X] = ode113(@eqs_motion, t_range, x_initial, options);

% PLotting
figure 
hold on
axis([-2 2 -0.5 0.5]);

planet_1 = plot(X(1, 1), X(1, 2), 'marker', '.', 'MarkerSize', 30, 'Color', 'b');
planet_2 = plot(X(1, 3), X(1, 4), 'marker', '.', 'MarkerSize', 30, 'Color', 'g');
planet_3 = plot(X(1, 5), X(1, 6), 'marker', '.', 'MarkerSize', 30, 'Color', 'r');
for i = 1:height(X)
    % Update positions of the planets
    set(planet_1, 'XData', X(i, 1), 'YData', X(i, 2));
    set(planet_2, 'XData', X(i, 3), 'YData', X(i, 4));
    set(planet_3, 'XData', X(i, 5), 'YData', X(i, 6));
    plot ( X(1:i, 1), X(1:i, 2), 'blue', ...
         X(1:i, 3), X(1:i, 4), 'green', ...
         X(1:i, 5), X(1:i, 6), 'red' )

    pause(0.0000001) 
end

% Legend for planets
title("Matlab's ode_113")
legend([planet_1, planet_2, planet_3], {'Planet 1', 'Planet 2', 'Planet 3'}, 'Location', 'southoutside');

%-------------------------------------------------------------------------------------------------------------------

% método RK4
h =0.0001; % Step size
n = (t_final-t_initial)/h; % Number of steps
steps = n-1;
t = linspace(t_initial, t_final, n);

u = [p1_start, p2_start, p3_start]; %Array with positions [x1,y1,x2,y2,x3,y3]
v = [v1_start, v2_start, v3_start]; %Array with velocities [vx1,vy1,vx2,vy2,x3,y3]
x = [u, v];

for i = 1:steps
    % Runge-Kutta 4th order method
    
    k1 = h .* v;
    k_1  = h .* f(u);

    k2 = h .* (v + k_1.*(1/2));
    k_2 = h .* f(u + k1/2);

    k3 = h .* (v + k_2/2);
    k_3 = h .* f(u + k2/2);

    k4 = h .* (v + k_3/2);
    k_4 = h .* f(u + k3/2);
    
    % Update state variable
    u = u + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
    v = v + (1/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
    x(i+1,:) = [u,v];

end

% Plotting
figure;
hold on
axis([-2 2 -0.5 0.5]);

planet_1RK4 = plot(x(1, 1), x(1, 2), 'marker', '.', 'MarkerSize', 30, 'Color', 'b');
planet_2RK4 = plot(x(1, 3), x(1, 4), 'marker', '.', 'MarkerSize', 30, 'Color', 'g');
planet_3RK4 = plot(x(1, 5), x(1, 6), 'marker', '.', 'MarkerSize', 30, 'Color', 'r');
for i = 1:100:height(x)
    % Update positions of the planets
    
    set(planet_1RK4, 'XData', x(i, 1), 'YData', x(i, 2));
    set(planet_2RK4, 'XData', x(i, 3), 'YData', x(i, 4));
    set(planet_3RK4, 'XData', x(i, 5), 'YData', x(i, 6));
    plot ( x(1:i, 1), x(1:i, 2), 'blue', ...
         x(1:i, 3), x(1:i, 4), 'green', ...
         x(1:i, 5), x(1:i, 6), 'red' )

    pause(0.0000001) 
end

% Legend for planets
title('Método RK4')
legend([planet_1RK4, planet_2RK4, planet_3RK4], {'planet_1', 'planet_2', 'planet_3'}, 'Location', 'southoutside');

%-------------------------------------------------------------------------------------------------------------------

%Método implícito de Euler con predictor-corrector

x1 = [u,v];
for i=1:steps
    %pred_x=x(i,:)+h.*x(i,:);
    pred_x=x1(i,:)+h*f_matrix(x1(i,:))';
    %u=x(i,:)+h.*pred_x;
    %v=x(i,:)+h.*f_matrix(pred_xt);
    %x(i+1,:)=[u,v];
    x1(i+1,:)=x1(i,:)+h*f_matrix(pred_x)';
end

% Plotting
figure;
hold on
axis([-2 2 -0.5 0.5]);
title('Metodo de Euler implicito con predictor-corrector')
planet_1E = plot(x1(1, 1), x1(1, 2), 'marker', '.', 'MarkerSize', 30, 'Color', 'b');
planet_2E = plot(x1(1, 3), x1(1, 4), 'marker', '.', 'MarkerSize', 30, 'Color', 'g');
planet_3E = plot(x1(1, 5), x1(1, 6), 'marker', '.', 'MarkerSize', 30, 'Color', 'r');
for i = 1:100:height(x1)
    % Update positions of the planets
    
    set(planet_1E, 'XData', x1(i, 1), 'YData', x1(i, 2));
    set(planet_2E, 'XData', x1(i, 3), 'YData', x1(i, 4));
    set(planet_3E, 'XData', x1(i, 5), 'YData', x1(i, 6));
    plot ( x1(1:i, 1), x1(1:i, 2), 'blue', ...
         x1(1:i, 3), x1(1:i, 4), 'green', ...
         x1(1:i, 5), x1(1:i, 6), 'red' )

    pause(0.0000001) 
end

% Legend for planets
title('Método de Euler implícito')
legend([planet_1E, planet_2E, planet_3E], {'planet_1', 'planet_2', 'planet_3'}, 'Location', 'southoutside');


%---------------------------------------------------------------------------------------------------------------
% RK4 method with matrix

h = 0.0001; % Step size
n = (t_final-t_initial)/h; % Number of steps
steps = n-1;
t = linspace(t_initial, t_final, n);

x = x_initial;

for i = 1:steps
    k1 = h * f_matrix(x(i, :));
    k2 = h * f_matrix(x(i, :) + k1'/2);
    k3 = h * f_matrix(x(i, :) + k2'/2);
    k4 = h * f_matrix(x(i, :) + k3');
    
    % Update state variable
    x(i+1, :) = x(i, :) + (1/6) * (k1 + 2*k2 + 2*k3 + k4)';
end

% Plotting
figure;
hold on
axis([-2 2 -0.5 0.5]);

planet_1RK4 = plot(x(1, 1), x(1, 2), 'marker', '.', 'MarkerSize', 30, 'Color', 'b');
planet_2RK4 = plot(x(1, 3), x(1, 4), 'marker', '.', 'MarkerSize', 30, 'Color', 'g');
planet_3RK4 = plot(x(1, 5), x(1, 6), 'marker', '.', 'MarkerSize', 30, 'Color', 'r');
for i = 1:100:height(x)
    % Update positions of the planets
    
    set(planet_1RK4, 'XData', x(i, 1), 'YData', x(i, 2));
    set(planet_2RK4, 'XData', x(i, 3), 'YData', x(i, 4));
    set(planet_3RK4, 'XData', x(i, 5), 'YData', x(i, 6));
    plot ( x(1:i, 1), x(1:i, 2), 'blue', ...
         x(1:i, 3), x(1:i, 4), 'green', ...
         x(1:i, 5), x(1:i, 6), 'red' )

    pause(0.0000001) 
end

% Legend for planets
title('Método RK4 con matrices')
legend([planet_1RK4, planet_2RK4, planet_3RK4], {'planet_1', 'planet_2', 'planet_3'}, 'Location', 'southoutside');

%------------------------------------------------------------------------------------------------------------------

function [xdot] = eqs_motion(~, x)
    % Unpack positions and velocities from x
    positions = x(1:6);
    velocities = x(7:12);
    
    % Find the accelerations acting on the planets
    accelerations = gravity_force(positions');
    
    % Pack accelerations into the output vector
    xdot = [velocities; accelerations'];
end

function y = f(x)
    positions = x(1:6);
    accelerations = gravity_force(positions);
    y = accelerations;
end

function y = f_matrix(x)
    %  We want to compute x_dot = f(t, x)  = Ax + b
    % Find b, which accounts for non linear contributions in the ode (in this case gravity force)
    positions = x(1:6);
    accelerations = gravity_force(positions);
    b = [zeros(1,6), accelerations];
    
    % Find A, which accounts for the linear part of ode 
    A = zeros(12);
    A(1:6, 7:end) = eye(6);
    
    % Compute f(t, x) = Ax + b
    y = A * x' + b';
end

% Define the accelerations function
function F = gravity_force(x)
    % define parameters
    global G
    global m1
    global m2
    global m3
    
    % Unpack positions for each planet (we assume that x=[x1,y1,x2,y2,x3,y3])
    p1 = x(1:2);
    p2 = x(3:4);
    p3 = x(5:6);

    % Compute accelerations for each planet
    a1 = -G * m2 * (p1 - p2) / norm(p1 - p2)^3 - G * m3 * (p1 - p3) / norm(p1 - p3)^3;
    a2 = -G * m3 * (p2 - p3) / norm(p2 - p3)^3 - G * m1 * (p2 - p1) / norm(p2 - p1)^3;
    a3 = -G * m1 * (p3 - p1) / norm(p3 - p1)^3 - G * m2 * (p3 - p2) / norm(p3 - p2)^3;

    F =[a1,a2,a3];
end
