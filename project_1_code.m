clear
clc
syms y
%% Constants 

% Wing Design Parameters
k = 0.25;
c1 = 0.35; c2 = 0.4; s = 2; % m 

% Wing Parameters
c = (c2+c1) + (((c1/2)-c2)/s)*y; % m
e = (((c2 + c1) + (((c1 / 2) - c2) / s) * y)*3/4)-c1; % m
A = s*((5/4)*c1 + (c2/2));
c_root = c1+c2;
c_tip = 1.5*c1;
c_MAC = double(((2/A)*int(c^2,0,s))); % m
c_mean = (2.5*c1+c2)/2;
eta = y/s; % 1/m 
C_L_a = 2*pi*sqrt(1-eta^2);
GJ = @(y) 8500*(1-(k*eta)); % Nm^2/rad
alpha_i = (5 - 3*eta)*(pi/180); % rad
n_values = []; % Store values of n for plotting
target_value = 0.1/100;

% Flight Conditions
U = 70; % m/s
rho = 1.225;    % Air density in kg/m^3
q_70 = 0.5 * rho * U^2; % Dynamic pressure at 70 m/s


%% Part A
% Elastic Wing Calculations for Divergance Dynamic Pressure 
% Wing Tip Deflection at 50 % q_d


% Initial value of n (Number of modes)
n = 1;

% Loop until eps reaches the target value
while true
    disp(n)

    % Create a cell array to store the vector of functions for current n
    functions = sym('f', [1 n]);
    diff_f = sym('f', [1 n]); 

    % Populate the vector with functions f_{i,j}(y) = j * y^i
    for i = 1:n
        functions(i) = i * y^i;
        diff_f(i) = diff(functions(i), y); % Store the derivative as a function vector 
    end

    % Initialize E matrix to store integrals
    E_matrix = sym('f', [n n]);
    for i = 1:n
        for j = 1:n 
            E_integrand = GJ*diff_f(i)*diff_f(j);
            E_matrix(i,j) = int(E_integrand, 0, s);  
        end
    end

    % Calculate K based on the functions vector
    K_matrix = sym('f', [n n]);
    for i = 1:n
        for j = 1:n
            K_integrand = (c^2)*e*C_L_a*functions(i)*functions(j);
            K_matrix(i,j) = int(-K_integrand, 0, s);  
        end
    end

    % Calculate F based on the functions vector
    F_column = sym('f', [n 1]);
    for i = 1:n
        F_integrand = (c^2)*e*C_L_a*alpha_i*functions(i);
        F_column(i) = int(F_integrand, 0, s);  
    end

        q_d = eig(double(E_matrix), -1*double(K_matrix));
        q_d = min(q_d(q_d >= 0));
        q_d_values(n) = q_d; % Store q_d for current n
        n_values(n) = n;     % Store the current mode number for plotting
     
    if n > 1   
        % Calculate error
        eps = abs((q_d - q_d_values(n-1)) / q_d_values(n-1));
        if eps <= target_value
            fprintf('\n_______SOLUTION FOUND_______')
            fprintf('\nTarget value reached with n = %d\n', n);
            fprintf('eps = %.5f\n', eps);
            fprintf('q_d = %.2f\n', q_d);
            break;
        end
    end
    n = n + 1;
end

% Plot the Variation of Dynamic Pressure with Number of Modes

% Interpolate 
n_dense = linspace(min(n_values), max(n_values), 100);  
q_d_dense = interp1(n_values, q_d_values, n_dense, 'pchip');  

figure;
plot(n_dense, q_d_dense, '-','LineWidth', 2);  
hold on;
plot(n_values, q_d_values, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r');  
xlabel('Number of Modes (n)');
ylabel('Dynamic Pressure q_d [N/m^2]');
title('Variation of Dynamic Pressure with Number of Modes');
grid on;


% wing tip deflection at 50% of divergence dynamic pressure

theta = inv(double(E_matrix) + .5*q_d*double(K_matrix)) * .5*q_d*double(F_column);
% for wing tip 
functions_tip = subs(functions,y,s);
for i = 1:n
    theta_t(i) = double(theta(i)*functions_tip(i));

end
theta_tip = sum(double(theta_t));
theta_tip_deg = theta_tip*(180/pi);

fprintf('\n Wing Tip Deflection at 50 percent q_d = %.2f degrees \n', theta_tip_deg);


%% Part B

% Guessed Wing Parameters
k = .25;
c1 = 0.24; c2 = 0.24; s = 1.28; % m 

% Updated Parameters
A = s*((5/4)*c1 + (c2/2));
c_root = c1+c2;
c_tip = 1.5*c1;
c_MAC = double(((2/A)*int(c^2,0,s))); % m
c_mean = (2.5*c1+c2)/2;
eta = y/s; % 1/m 
C_L_a = 2*pi*sqrt(1-eta^2);
GJ = @(y) 8500*(1-(k*eta)); % Nm^2/rad
alpha_i = (5 - 3*eta)*(pi/180); % rad

n_values = []; 
target_value = 0.1/100;

% Main Loop
% Initial value of n (Number of modes)
n = 6;

% Loop until eps reaches the target value
while true

    % Create a cell array to store the vector of functions for current n
    functions = sym('f', [1 n]);
    diff_f = sym('f', [1 n]); 

    % Populate the vector with functions f_{i,j}(y) = j * y^i
    for i = 1:n
        functions(i) = i * y^i;
        diff_f(i) = diff(functions(i), y); % Store the derivative as a function vector 
    end

    % Initialize E matrix to store integrals
    E_matrix = sym('f', [n n]);
    for i = 1:n
        for j = 1:n 
            E_integrand = GJ*diff_f(i)*diff_f(j);
            E_matrix(i,j) = int(E_integrand, 0, s);  
        end
    end

    % Calculate K based on the functions vector
    K_matrix = sym('f', [n n]);
    for i = 1:n
        for j = 1:n
            K_integrand = (c^2)*e*C_L_a*functions(i)*functions(j);
            K_matrix(i,j) = int(-K_integrand, 0, s);  
        end
    end
    
    % Calculate F based on the functions vector
    F_column = sym('f', [n 1]);
    for i = 1:n
        F_integrand = (c^2)*e*C_L_a*alpha_i*functions(i);
        F_column(i) = int(F_integrand, 0, s);  
    end

    if n > 1
        q_d = eig(double(E_matrix), -1*double(K_matrix));
        q_d = min(q_d(q_d>=0));
        q_d_values(n) = q_d; % Store q_d for current n

        % Calculate error
        eps = abs((q_d - q_d_values(n-1))/q_d_values(n-1));
        if eps <= target_value
            
            break;
        end
    end
    n = n + 1;
end

% Calculate theta at 70 m/s: 
theta = inv(double(E_matrix) + q_70*double(K_matrix)) * q_70*double(F_column);
% for wing tip 
functions_tip = subs(functions,y,s);
for i = 1:n
    theta_t(i) = double(theta(i)*functions_tip(i));

end
theta_tip = sum(double(theta_t));
theta_tip_deg = theta_tip*(180/pi);

% Calculate theta(y) and lift distribution along the span

y_values = linspace(0, s, 100);  
theta_y = zeros(size(y_values)); 
L_y = zeros(size(y_values));     
alpha_i_y = zeros(size(y_values));  

% Calculate theta(y) and lift distribution along the span
for j = 1:length(y_values)
    % Calculate local twist theta(y)
    functions_y = subs(functions, y, y_values(j));
    theta_y_temp = 0;
    for i = 1:n
        theta_y_temp = theta_y_temp + double(theta(i) * functions_y(i));
    end
    theta_y(j) = theta_y_temp;
    
    % Calculate local lift L'(y)
    C_L_y = double(subs(C_L_a * (alpha_i + theta_y_temp), y, y_values(j))); % Local lift coefficient
    L_y(j) = q_70 * A * C_L_y; % Local lift per unit span

    % Local chord
    c_y(j) = double(subs(c, y, y_values(j)));          

end


% Convert theta from radians to degrees
theta_y_deg = theta_y * (180/pi);

% Wing Tip Deflection
functions_tip = subs(functions, y, s);
for i = 1:n
    theta_t(i) = double(theta(i)*functions_tip(i));
end
theta_tip = sum(double(theta_t));
theta_tip_deg = theta_tip*(180/pi);

% Divergence 
U_d = sqrt((2*q_d)/rho);

% Total Lift and Moment (Elastic Wing)

L_total = double(trapz(y_values, L_y)); 
M_total = L_total*(s/3); % Elliptical wing distribution

% Output

fprintf('\nk = %.2fm, c1 = %.2fm  c2 = %.2fm and s = %.2fm\n', k, c1, c2, s);
fprintf('-------------------------------------------------\n');

if abs(theta_tip_deg) <= 1
    fprintf('Wing Tip Deflection = %.2f\n', theta_tip_deg);
else
    fprintf('Wing Tip Deflection is out of acceptable range (greater than 1 degree).\n');
end

if abs(U_d) >= 150
    fprintf('Divergence Speed = %.2f\n', U_d);
else
    fprintf('Divergence Speed is out of acceptable range (less than 150 m/s).\n');
end

if abs(L_total) >= 700
    fprintf('Lift = %.2f\n', L_total);
else
    fprintf('Lift is out of acceptable range (less than 700 N).\n');
    fprintf('Lift = %.2f\n', L_total);
end

if abs(M_total) <= 300
    fprintf('Moment = %.2f\n', M_total);
else
    fprintf('Moment is out of acceptable range (greater than 300 Nm).\n');
    fprintf('Moment = %.2f\n', M_total);
end 
    
if abs(s/c_mean) >= 3
    fprintf('s/c_mean = %.2f\n', s/c_mean);
else
    fprintf('s/c_mean is out of acceptable range (less than 3).\n');
    fprintf('s/c_mean = %.2f\n', s/c_mean);
end 
  
if abs(theta_tip_deg) <= 1 && abs(U_d) >= 150 && abs(L_total) >= 700 && abs(M_total) <= 300 && abs(s/c_mean) >= 3
    fprintf('----------------------------------\n');
    fprintf('_____All values are acceptable____\n');
end


%% Part C (uses Part B)

% Plot twist variation
figure;
plot(y_values, theta_y_deg, 'LineWidth', 2);
xlabel('Spanwise Location (y) [m]');
ylabel('Twist Angle \theta (degrees)');
title('Variation of Twist Angle Along the Span');
grid on;

%% Part D (uses Part B)
% elastic and rigid lift distribution along the span

% Rigid Wing Lift Distribution (No Twist)

L_y_rigid = zeros(size(y_values)); 

% Calculate lift distribution for a rigid wing
for j = 1:length(y_values)
    
    C_L_rigid = double(subs(C_L_a * alpha_i, y, y_values(j))); 
    L_y_rigid(j) = q_70 * A* C_L_rigid;            
end


% Plot 
plot(y_values, L_y, 'LineWidth', 2);  % Elastic Wing
hold on;
plot(y_values, L_y_rigid, 'LineWidth', 2);  % Rigid Wing

xlabel('Spanwise Location (y) [m]');
ylabel('Lift per unit span L''(y) [N/m]');
title('Lift Distribution for Elastic vs Rigid Wing');
legend('Elastic Wing', 'Rigid Wing');

grid on;
hold off;


%% Part E (uses Part B) 
% percentage increase in total lift and wing root bending moment


% Total Lift and Moment (Rigid Wing)

L_total_rigid = double(trapz(y_values, L_y_rigid)); % Numerical integration of lift along span
M_total_rigid = double((s/3)* L_total_rigid);


% Percentage increase in total lift and wing root bending moment 

P_L_inc = ((L_total - L_total_rigid)/L_total)*100;
P_M_inc = ((M_total - M_total_rigid)/M_total)*100;

% Output Results

fprintf('\nPercentage Increase in Total Lift and Moment = %.2f percent \n', P_L_inc);



%% Part F (uses Part B) 
% aileron reversal dynamic pressure and the rolling power of the wing

% Flight Conditions
U = 70; % m/s
rho = 1.225;    % Air density in kg/m^3
q_70 = 0.5 * rho * U^2; % Dynamic pressure at 70 m/s

% Aileron 
y1 = 0.5*s; y2 = 0.8*s; 

% Aileron chord
c_a = .3*c_y(65);

phi = acos(2*c_a-1);
C_L_beta = ((-sin(phi))*(1-cos(phi)))/2;
C_M_beta = 2*(pi-phi+sin(phi));

f = y*(4*s-2*y);
f_d = diff(f);

E_integrand_F = GJ*f_d*f_d;
E = int(E_integrand_F,0,s);

C_integrand_F = c*C_L_a*f*y;
C = int(C_integrand_F,0,s);

A_integrand_F = c_a*C_L_beta*y;
A = int(A_integrand_F,y1,y2);

D_integrand_F = (c^2)*e*C_L_a*f*f;
D = int(D_integrand_F,0,s);

B_integrand_F = (c_a^2)*((e*C_L_beta)+C_M_beta)*f;
B = int(B_integrand_F,y1,y2);

F_integrand_F = c_a*C_L_a*y^2;
F = int(F_integrand_F,0,s);

G_integrand_F = (c_a^2)*e*C_L_a*f*y;
G = int(G_integrand_F,0,s);

q_r = double((A*E)/((A*D)-(B*C)));
R_p = -1*double(U*   ((A*(E/q_70-D))+(B*C))  /   (  (C*G)  +  (F*  ((E/q_70)-D) )))  ;


% Output Results

fprintf('\nAileron Reversal Dynamic Pressure = %.2f \n', q_r);
fprintf('Rolling Power at 70 m/s = %.2f \n', R_p);

%% Part G: Study effects of Chord Size on Rolling Power

% Wing Parameters
k = .25;
c1 = 0.35; c2 = 0.4; s = 2; % m 

% Parameters
c = (c2+c1) + (((c1/2)-c2)/s)*y; % m
e = (((c2 + c1) + (((c1 / 2) - c2) / s) * y)*3/4)-c1; % m
A = s*((5/4)*c1 + (c2/2));
c_root = c1+c2;
c_tip = 1.5*c1;
c_MAC = double(((2/A)*int(c^2,0,s))); % m
c_mean = (2.5*c1+c2)/2;
eta = y/s; % 1/m 
C_L_a = 2*pi*sqrt(1-eta^2);
GJ = @(y) 8500*(1-(k*eta)); % Nm^2/rad
alpha_i = (5 - 3*eta)*(pi/180); % rad
% Flight Conditions
U = 70; % m/s
rho = 1.225;    
q_70 = 0.5 * rho * U^2; 

% Aileron 
y1 = 0.5*s; y2 = 0.8*s; 

% Aileron at different chords
c_a_values = [.3, .4, .5, .6];

R_p_results = zeros(1, length(c_a_values));

for i = 1:length(c_a_values)

    c_a = c_a_values(i)*c_y(65);

    phi = acos(2*c_a-1);

    C_L_beta = ((-sin(phi))*(1-cos(phi)))/2;
    C_M_beta = 2*(pi-phi+sin(phi));

    f = y*(4*s-2*y);
    f_d = diff(f);

    E_integrand_F = GJ*f_d*f_d;
    E = int(E_integrand_F,0,s);

    C_integrand_F = c*C_L_a*f*y;
    C = int(C_integrand_F,0,s);

    A_integrand_F = c_a*C_L_beta*y;
    A = int(A_integrand_F,y1,y2);

    D_integrand_F = (c^2)*e*C_L_a*f*f;
    D = int(D_integrand_F,0,s);

    B_integrand_F = (c_a^2)*((e*C_L_beta)+C_M_beta)*f;
    B = int(B_integrand_F,y1,y2);

    F_integrand_F = c_a*C_L_a*y^2;
    F = int(F_integrand_F,0,s);

    G_integrand_F = (c_a^2)*e*C_L_a*f*y;
    G = int(G_integrand_F,0,s);

    q_r = double((A*E)/((A*D)-(B*C)));
    R_p = 1*double(U*   ((A*(E/q_70-D))+(B*C))  /   (  (C*G)  +  (F*  ((E/q_70)-D) )))  ;

    R_p_results(i) = R_p;

    % Output Results
    fprintf('Rolling Power (Aileron Chord = %.1fc) = %.2f \n',c_a_values(i), R_p);
    

end

fprintf('---------------------------------------------\n');

% Plot the rolling power results
figure;
plot(c_a_values, R_p_results, '-o', 'LineWidth', 2);
xlabel('Aileron Chord Ratio (c_a / c)');
ylabel('Rolling Power (R_p) [rad/degree*s]');
title('Rolling Power vs Aileron Chord');
grid on;
    




