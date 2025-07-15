clear
clc
syms U R R_new R_new_m R_new_d R_best
%%
%%%%%%%%% Given Parameters %%%%%%%%%
c = 0.5;        % meters
k_1 = 5000;     % kN/m
k_theta1 = 500; % N.m/rad
k_2 = 1000;
c_theta1 = 0;
c_2 = 0;
m_wing = 5;
I_CG = 0.05;
x_g = 0.15;
m_1 = 2;
% kN/m
% N.m.s/rad
% N.s/rad
% kg
% kg.m^2
% meters
% kg
Appendix B
 x_m = 0.15;
U_max = 100;
rho = 1.225;
b = c/2;
% meters
% m/s
% kg/m^3
% meters
% meters
% distance fom leading edge to spring 1 in
s = 1;
x_1 = 0;
meters
x_2 = (3/4)*c; % distance fom leading edge to spring 2 in meters
%%%%%%%% Matrices %%%%%%%%%%
% Mass Matrix
m_total = m_wing + m_1;
S = -m_1*(b - x_m) + m_wing*(x_g - b);
I_E = I_CG + m_1*(b - x_m)^2 + m_wing*(x_g - b)^2;
M = [m_total S ; S I_E];
% Stiffness Matrix
k_h = k_1 + k_2;
k_theta = k_1*(b-x_1)^2 + k_2*(x_2 - b)^2 + k_theta1;
E_k = -k_1*b + k_2*(b/2);
E = [k_h E_k ; E_k k_theta];
% K Matrix
K = [0 1 ; 0 -(b/2)];
% Structural Damping Matrix
B_s = [c_2 c_2*(b/2) ; c_2*(b/2) c_2*(b/2)^2+c_theta1];
% Aerodynamic Damping Matrix
B_a = [1 (b/2) ; -(b/2) 0];
% B bar
B_bar = (pi*rho*s*c*U)*(B_a)+B_s;
% K bar
K_bar = (pi*rho*s*c*U^2)*K;
%characteristic eqn
ceqn = (R^2)*M + R*B_bar + K_bar + E;
det_ceqn = det(ceqn);
coef = coeffs(det_ceqn,R);
P_0 = coef(:,1);
P_1 = coef(:,2);
P_2 = coef(:,3);
P_3 = coef(:,4);

 P_4 = coef(:,5);
%% Q2 - Loop to Find Flutter Speed
eq = R^4*P_4 + R^3*P_3 + R^2*P_2 + R*P_1 + P_0 ;
U_values = 1:1:100;
V_critical = NaN;
found = false;
for i=1:length(U_values)
   U_current = U_values(i);
   eq_sub = subs(eq, U, U_current);
   roots_eq = double(solve(eq_sub, R));
   if any(real(roots_eq) > 0)
       V_critical = U_current;
       found = true;
break;
end end
% Output Flutter Speed
if found
fprintf ('Critical Velocity where the real part of the
root becomes positive is %.2f m/s \n',V_critical)
else
fprintf ('No positive real part detected up to Maximum Velocity %.2f m/s \n ', U_max)
end
%% Q3 - Plot Eigenvalues vs Airspeed
U_critical_new = V_critical*1.2;
U_values_2 = linspace(0, U_critical_new, 500); % Define range of U
real_parts_1 = cell(1, length(U_values_2));
imag_parts_1 = cell(1, length(U_values_2));
real_parts_2 = cell(1, length(U_values_2));
imag_parts_2 = cell(1, length(U_values_2));
real_parts_3 = cell(1, length(U_values_2));

 imag_parts_3 = cell(1, length(U_values_2));
real_parts_4 = cell(1, length(U_values_2));
imag_parts_4 = cell(1, length(U_values_2));
for i = 1:length(U_values_2)
   U_current_2 = U_values_2(i);
   eq_sub = subs(eq, U, U_current_2);
   eigenvalues = double(solve(eq_sub, R));
   real_parts_1{i} = real(eigenvalues(1));  % Store all real
parts
imag_parts_1{i} = imag(eigenvalues(1)); % Store all imaginary parts
   real_parts_2{i} = real(eigenvalues(2));  % Store all real
parts
   imag_parts_2{i} = imag(eigenvalues(2));
   real_parts_3{i} = real(eigenvalues(3));  % Store all real
parts
   imag_parts_3{i} = imag(eigenvalues(3));
   real_parts_4{i} = real(eigenvalues(4));  % Store all real
parts
   imag_parts_4{i} = imag(eigenvalues(4));
end
% Combine into matrices for plotting (one row for each U value)
real_parts_matrix_1 = cell2mat(cellfun(@(x) real_parts_1, 'UniformOutput', false)); imag_parts_matrix_1 = cell2mat(cellfun(@(x) imag_parts_1, 'UniformOutput', false)); real_parts_matrix_2 = cell2mat(cellfun(@(x) real_parts_2, 'UniformOutput', false)); imag_parts_matrix_2 = cell2mat(cellfun(@(x) imag_parts_2, 'UniformOutput', false)); real_parts_matrix_3 = cell2mat(cellfun(@(x) real_parts_3, 'UniformOutput', false));
x(:)',
x(:)',
x(:)',
x(:)',
x(:)',

 imag_parts_matrix_3 = cell2mat(cellfun(@(x) imag_parts_3, 'UniformOutput', false)); real_parts_matrix_4 = cell2mat(cellfun(@(x) real_parts_4, 'UniformOutput', false)); imag_parts_matrix_4 = cell2mat(cellfun(@(x) imag_parts_4, 'UniformOutput', false));
x(:)',
x(:)',
x(:)',
% Prepare Airspeed values for each eigenvalue (same size as real_parts_matrix)
U_matrix = repmat(U_values_2(:), 1, size(real_parts_matrix_1, 2));
% Plot Real Parts vs Airspeed (U_values_2)
figure;
plot(U_matrix, real_parts_matrix_1, 'r.', 'MarkerSize', 10); hold on % Plot real parts
plot(U_matrix, real_parts_matrix_2, 'g.', 'MarkerSize', 10); plot(U_matrix, real_parts_matrix_3, 'b.', 'MarkerSize', 10); plot(U_matrix, real_parts_matrix_4, 'y.', 'MarkerSize', 10); xline(V_critical, '--k', 'Critical Speed', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Airspeed (U) [m/s]');
ylabel('Real Part of Eigenvalues');
title('Real Part of Eigenvalues vs Airspeed');
% Plot Imaginary Parts vs Airspeed (U_values_2)
figure;
plot(U_matrix, imag_parts_matrix_1, 'r.', 'MarkerSize', 10); hold on % Plot imaginary parts
plot(U_matrix, imag_parts_matrix_2, 'g.', 'MarkerSize', 10); plot(U_matrix, imag_parts_matrix_3, 'b.', 'MarkerSize', 10); plot(U_matrix, imag_parts_matrix_4, 'y.', 'MarkerSize', 10); xline(V_critical, '--k', 'Critical Speed', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Airspeed (U) [m/s]');
ylabel('Imaginary Part of Eigenvalues');
title('Imaginary Part of Eigenvalues vs Airspeed');
%% Q4 - Variation of Critical Speed vs Spring Constants (Hard Coded)

 % 5th iteration
k_1 = 5000*5;     % kN/m
k_theta1 = 500*5; % N.m/rad
k_2 = 1000*5;     % kN/m
%%%%%%%% Matrices %%%%%%%%%%
% Mass Matrix
m_total = m_wing + m_1;
S = -m_1*(b - x_m) + m_wing*(x_g - b);
I_E = I_CG + m_1*(b - x_m)^2 + m_wing*(x_g - b)^2;
M = [m_total S ; S I_E];
% Stiffness Matrix
k_h = k_1 + k_2;
k_theta = k_1*(b-x_1)^2 + k_2*(x_2 - b)^2 + k_theta1;
E_k = -k_1*b + k_2*(b/2);
E = [k_h E_k ; E_k k_theta];
% K Matrix
K = [0 1 ; 0 -(b/2)];
% Structural Damping Matrix
B_s = [c_2 c_2*(b/2) ; c_2*(b/2) c_2*(b/2)^2+c_theta1];
% Aerodynamic Damping Matrix
B_a = [1 (b/2) ; -(b/2) 0];
% B bar
B_bar = (pi*rho*s*c*U)*(B_a)+B_s;
% K bar
K_bar = (pi*rho*s*c*U^2)*K;
%characteristic eqn
ceqn = (R^2)*M + R*B_bar + K_bar + E;
det_ceqn = det(ceqn);
coef = coeffs(det_ceqn,R);
P_0 = coef(:,1);
P_1 = coef(:,2);
P_2 = coef(:,3);
P_3 = coef(:,4);
P_4 = coef(:,5);
eq = R^4*P_4 + R^3*P_3 + R^2*P_2 + R*P_1 + P_0 ;
U_values = 1:1:100;
V_critical = NaN;
found = false;

 for i=1:length(U_values)
   U_current = U_values(i);
   eq_sub = subs(eq, U, U_current);
   roots_eq = double(solve(eq_sub, R));
   if any(real(roots_eq) > 0)
       V_critical = U_current;
       found = true;
break;
end end
% Output Flutter Speed
if found
fprintf ('Critical Velocity where the real part of the
root becomes positive is %.2f m/s \n',V_critical_5)
else
fprintf ('No positive real part detected up to Maximum Velocity %.2f m/s \n ', U_max)
end
%% Q4 - Plot
V_critical_v = [0; V_critical; V_critical_2; V_critical_3; V_critical_4; V_critical_5]
magnitudes = [0; 1; 2; 3; 4; 5]
figure;
plot(magnitudes, V_critical_v, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Magnitude of K''s');
ylabel('Critical Velocity (m/s)');
title('Critical Velocity vs Magnitude of K''s');
%% Q5 - Best Combination of K_1 K_2 and K_theta1 for Maximum
Critical Speed and minimum stiffness
% Create ranges for k1, k2, and k_theta1
k1_range_best = 1000:1:9000;          % 1 kN/m to 9 kN/m

 k2_range_best = 1000:1:9000;          % 1 kN/m to 9 kN/m
k_theta1_range_best = 0:1:700;         % 0 to 700 Nm/rad
% Initialize variables for tracking results
max_critical_speed = 0;
best_combination = [0, 0, 0];
critical_speeds_best = zeros(length(k1_range_best), 1); %%%%%%%%% remove
% Iterate through all combinations of k1, k2, and k_theta1 for i = 1:length(k1_range_best)
   k1_best = k1_range_best(i);
   for j = 1:length(k2_range_best)
       k2_best = k2_range_best(j);
if k1_best + k2_best <= 10000 % Enforce sum of spring constants constraint
           for k = 1:length(k_theta1_range_best)
               k_theta1_best = k_theta1_range_best(k);
                   if k_theta1_best <= 700
                       %Change Affected Matrix E
k_h_best = k1_best + k2_best; k_theta_best = k1_best*(b-x_1)^2 +
k2_best*(x_2 - b)^2 + k_theta1_best;
E_k_best = -k1_best*b +
k2_best*(b/2);
E_k_best k_theta_best];
E_best = [k_h_best E_k_best ;
R_best*B_bar + K_bar + E_best;
% New characteristic equation
ceqn_best = (R_best^2)*M + det_ceqn_best = det(ceqn_best);

 R_best);
coef_best = coeffs(det_ceqn_best,
if length(coef_best) <5
coef_best = [zeros(1,
5-length(coef_best)), coef_best];
end
                       P_0_best = coef_best(:,1);
                       P_1_best = coef_best(:,2);
                       P_2_best = coef_best(:,3);
                       P_3_best = coef_best(:,4);
                       P_4_best = coef_best(:,5);
C_eq_best = R_best^4*P_4_best + R_best^3*P_3_best + R_best^2*P_2_best + R_best*P_1_best +
P_0_best ;
U_values_best(l);
U_current_best);
double(solve(eq_sub_best, R_best));
U_current_best;
if any(real(roots_eq_best) > 0) V_critical_best =
found = true;
break; end
U_values_best = 0:1:50;
found = false;
V_critical_best = NaN;
for l=1:length(U_values_best)
U_current_best = eq_sub_best = subs(C_eq_best, U, roots_eq_best =

 end end
               % Update if a higher critical speed is found
                if ~isnan(V_critical_best) & V_critical_best
> max_critical_speed
max_critical_speed = V_critical_best; best_combination = [k1_best, k2_best,
k_theta1_best];
end end
end end
end
% Output results
fprintf('Maximum
max_critical_speed);
fprintf('Optimal Spring Constants:\n');
fprintf('  k_1: %.2f kN/m\n', best_combination(1) / 1000);
fprintf('  k_2: %.2f kN/m\n', best_combination(2) / 1000);
fprintf('  k_{theta1}: %.2f Nm/rad\n', best_combination(3));
%% Q6 - Critical Speed vs x_m
% Range for x_m (point mass location)
xm_range = 0:.01:c;
% Initialize results
critical_speeds_m = zeros(length(xm_range), 1);
for i = 1:length(xm_range)
   x_m_current = xm_range(i); % Current point mass location
   % New Mass Matrix
   m_total_new = m_wing + m_1;
   S_new = -m_1*(b - x_m_current) + m_wing*(x_g - b);
I_E_new = I_CG + m_1*(b - x_m_current)^2 + m_wing*(x_g - b)^2;
   M_new = [m_total_new S_new ; S_new I_E_new];
%%%%%%%%%
Critical
Speed: %.2f
m/s\n',

    % Characteristic equation for New Mass
ceqn_new_m = (R_new_m^2)*M_new + R_new_m*B_bar + K_bar + E;
   det_ceqn_new_m = det(ceqn_new_m);
   coef_new_m = coeffs(det_ceqn_new_m, R_new_m);
if length(coef_new_m) <5
coef_new_m = [zeros(1, 5-length(coef_new_m)),
coef_new_m];
end
   P_0_m = coef_new_m(:,1);
   P_1_m = coef_new_m(:,2);
   P_2_m = coef_new_m(:,3);
   P_3_m = coef_new_m(:,4);
   P_4_m = coef_new_m(:,5);
C_eq_new_m = R_new_m^4*P_4_m + R_new_m^3*P_3_m + R_new_m^2*P_2_m + R_new_m*P_1_m + P_0_m ;
   U_values_m = 0:1:500;
   found = false;
   for j=1:length(U_values_m)
       U_current_m = U_values_m(j);
       eq_sub_m = subs(C_eq_new_m, U, U_current_m);
       roots_eq_m = double(solve(eq_sub_m, R_new_m));
       if any(real(roots_eq_m) > 0)
        critical_speeds_m(i) = U_current_m;
           found = true;
       break;
end end
end
% Plot the results

 figure;
plot(xm_range, critical_speeds_m, 'r-', 'LineWidth', 2); xlabel('Point Mass Location x_m (m)');
ylabel('Critical Speed (m/s)');
title('Variation of Critical Speeds vs Point Mass Location');
grid on;
%%%%%%% Find the optimal location for maximum critical speed %%%%%%%%%%%%%
[max_speed, max_idx] = max(critical_speeds_m);
optimal_xm = xm_range(max_idx);
fprintf('Maximum Critical Speed: %.2f m/s\n', max_speed); fprintf('Optimal Point Mass Location (x_m): %.2f m\n', optimal_xm);
%% Q7 - Variation of Critical Speed Versus the Damping Constants
% Define ranges for c2 and c_theta1
c_2_range = linspace(0, 20, 20);
c_theta1_range = linspace(0, 0.7, 20);
% Initialize results
critical_speeds_c2 = zeros(length(c_2_range), 1); critical_speeds_c_theta1 = zeros(length(c_theta1_range), 1); % Symbolic variables
syms R_new_d U
% Fix constants for each loop
fixed_c_theta1 = c_theta1_range(1); % Constant c_theta1 when varying c_2
fixed_c_2 = c_2_range(1); % Constant c_2 when varying c_theta1
% Loop over c2_range (varying c_2 while keeping c_theta1 constant)
for i = 1:length(c_2_range)
   c_2_new = c_2_range(i);
   c_theta1_new = fixed_c_theta1; % Fix c_theta1
   % Update affected matrix B_s
B_s_new = [c_2_new c_2_new*(b/2); c_2_new*(b/2) c_2_new*(b/2)^2 + c_theta1_new];
   B_bar_new = (pi*rho*s*c*U)*B_a + B_s_new;

    % Characteristic equation for new damping
ceqn_new_d = (R_new_d^2)*M + R_new_d*B_bar_new + K_bar + E;
   det_ceqn_new_d = det(ceqn_new_d);
   coef_new_d = coeffs(det_ceqn_new_d, R_new_d);
   % Ensure coefficients are 5 terms
coef_new_d = [zeros(1, 5-length(coef_new_d)), coef_new_d];
   % Construct characteristic polynomial
C_eq_new_d = coef_new_d(5)*R_new_d^4 + coef_new_d(4)*R_new_d^3 + ...
                                  coef_new_d(3)*R_new_d^2  +
coef_new_d(2)*R_new_d + coef_new_d(1);
   % Iterate over possible U values
   U_values_d = linspace(0, 500, 500);
   for j = 1:length(U_values_d)
       U_current_d = U_values_d(j);
       eq_sub_d = subs(C_eq_new_d, U, U_current_d);
       roots_eq_d = double(solve(eq_sub_d, R_new_d));
       % Find critical speed
       if any(real(roots_eq_d) > 0)
           critical_speeds_c2(i) = U_current_d;
break; end
end end
% Loop over c_theta1_range (varying c_theta1 while keeping c_2 constant)
for k = 1:length(c_theta1_range)
   c_2_new = fixed_c_2; % Fix c_2
   c_theta1_new = c_theta1_range(k);
   % Update affected matrix B_s
B_s_new = [c_2_new c_2_new*(b/2); c_2_new*(b/2) c_2_new*(b/2)^2 + c_theta1_new];
   B_bar_new = (pi*rho*s*c*U)*B_a + B_s_new;
   % Characteristic equation for new damping
ceqn_new_d = (R_new_d^2)*M + R_new_d*B_bar_new + K_bar + E;

    det_ceqn_new_d = det(ceqn_new_d);
   coef_new_d = coeffs(det_ceqn_new_d, R_new_d);
   % Ensure coefficients are 5 terms
coef_new_d = [zeros(1, 5-length(coef_new_d)), coef_new_d];
   % Construct characteristic polynomial
C_eq_new_d = coef_new_d(5)*R_new_d^4 + coef_new_d(4)*R_new_d^3 + ...
                                  coef_new_d(3)*R_new_d^2  +
coef_new_d(2)*R_new_d + coef_new_d(1);
   % Iterate over possible U values
   U_values_d = linspace(0, 500, 500);
   for j = 1:length(U_values_d)
       U_current_d = U_values_d(j);
       eq_sub_d = subs(C_eq_new_d, U, U_current_d);
       roots_eq_d = double(solve(eq_sub_d, R_new_d));
       % Find critical speed
       if any(real(roots_eq_d) > 0)
           critical_speeds_c_theta1(k) = U_current_d;
break; end
end end
% Plot critical speeds vs c_2
figure;
plot(c_2_range, critical_speeds_c2, 'b-', 'LineWidth', 1.5); xlabel('Parameter Variation c_2');
ylabel('Critical Speed (m/s)');
grid on;
title(sprintf('Critical Speeds vs c_2
%.2f)', fixed_c_theta1));
% Plot critical speeds vs c_theta1
figure;
plot(c_theta1_range, critical_speeds_c_theta1, 'r-', 'LineWidth', 1.5);
xlabel('Parameter Variation c_{\theta1}');
ylabel('Critical Speed (m/s)');
grid on;
(c_{\\theta1} =

title(sprintf('Critical Speeds vs c_{\\theta1} (c_2 = %.2f)', fixed_c_2));
