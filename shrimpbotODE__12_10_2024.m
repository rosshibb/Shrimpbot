% Constants

%physical leg constants:
M = 0.0183;  % Leg assembly mass (kg)
I_c = 18106.53e-9; %moment of intertia of leg assembly about an axis running through the center of mass and parallel to the leg's pivot axis (kg*m^2)
l_c = 0.03897;  % Distance from leg COM to leg pivot axis (m)
rho = 997;  % Water density (kg/m^3)
C_D_r = 1;  % Drag coefficient for the rigid part of the leg
% C_D_c = 0; %Drag coefficient for the compliant part of the leg; will change depending on direction of motion
%distances along the leg from pivot axis to ends of each of the 6 regions of the leg (m):
h1 = 0.0093;
h2 = 0.0293;
h3 = 0.0415;
h4 = 0.0583;
h5 = 0.0755;
h6 = 0.0841;
%widths of the ends of each of the 6 regions of the leg (m):
w1 = 0.0186;
w2 = 0.0093;
w3 = 0.0233;
w4 = 0.0360;
w5 = 0.0388;
w6 = 0.0141;
%solenoid constants:
m = 0.6484;  % Magnetic dipole moment of permanent magnet stack (A*m^2)
N = 300;  % Number of turns of coil in solenoid
mu = 4 * pi * 1e-7;  % Absolute electromagnetic permeability of water (H/m)
L = 0.0155; % Length of solenoid (m)
%constants for the shrimp metachronal swimming theta function
a = 40;
b = -20;
c = 0.55;
d_P_S = 0.8;
d_R_S = 1.1;


% Functions for theta(t) (piecewise, depending on t) to model shrimp 
% metachronal swimming, along with their first and second derivatives:
theta1 = @(t) b + d_P_S*(a/2) - (a/2)*cos(t*(pi/c)) - (a/2)*(d_P_S-1)*cos(t*((2*pi)/c));
theta1_dot = @(t) (a/2)*(pi/c)*sin(t*(pi/c)) + (a/2)*(d_P_S-1)*((2*pi)/c)*sin(t*((2*pi)/c));
theta1_dot_dot = @(t) (a/2)*(pi/c)^2*cos(t*(pi/c)) + (a/2)*(d_P_S-1)*((2*pi)/c)^2*cos(t*((2*pi)/c));

theta2 = @(t) b + d_R_S*(a/2) - (a/2)*cos((t-1)*(pi/(c-1))) - (a/2)*(d_R_S-1)*cos((t-1)*((2*pi)/(c-1)));
theta2_dot = @(t) (a/2)*((pi/(c-1)))*sin((t-1)*(pi/(c-1))) + (a/2)*(d_R_S-1)*(((2*pi)/(c-1)))*sin((t-1)*((2*pi)/(c-1)));
theta2_dot_dot = @(t) (a/2)*((pi/(c-1)))^2*cos((t-1)*(pi/(c-1))) + (a/2)*(d_R_S-1)*(((2*pi)/(c-1)))^2*cos((t-1)*((2*pi)/(c-1)));

% Functions isolating the current [I(t)] term in the ODE (there are two in
% order to account for the piecewise nature of the theta function
% I_iso = @(t) (((1/12)*M*l^2 + M*l_c^2) * theta(t) + ...
%               (1/8)*rho*C_d*(h^4*(w2-w1) + w1*l^4) * theta_dot(t)^2) ...
%               / (m*mu*n*L * cosd(theta(t)));
I_iso1 = @(t,C_D_c) ( (I_c + M*l_c^2) * (theta1_dot_dot(t))*(pi/180) + ...
            sign(theta1_dot(t))*(pi/180)^2*0.5*rho*theta1_dot(t)^2*(...
            (C_D_r*w1*((h1^4)/2))+...                                                   Region 1
            (C_D_r*w2*((h2^4-h1^4)/4))+...                                              Region 2
            ((C_D_r*(((w3-w2)/(h3-h2))*((4*h3^5-5*h2*h3^4+h2^5)/(20))+(w2*(h3-h2)))))+... Region 3
            ((C_D_c*(((w4-w3)/(h4-h3))*((4*h4^5-5*h3*h4^4+h3^5)/(20))+(w3*(h4-h3)))))+... Region 4
            ((C_D_c*(((w5-w4)/(h5-h4))*((4*h5^5-5*h4*h5^4+h4^5)/(20))+(w4*(h5-h4)))))+... Region 5
            ((C_D_c*(((w6-w5)/(h6-h5))*((4*h6^5-5*h5*h6^4+h5^5)/(20))+(w5*(h6-h5)))))...  Region 6
            )...
            )... 
            / (m*mu*N*(1/L) * cosd(theta1(t)));
 I_iso2 = @(t,C_D_c) ( (I_c + M*l_c^2) * theta2_dot_dot(t)*(pi/180) + ...
            sign(theta2_dot(t))*(pi/180)^2*0.5*rho*theta2_dot(t)^2*(...
            (C_D_r*w1*((h1^4)/2))+...                                                   Region 1
            (C_D_r*w2*((h2^4-h1^4)/4))+...                                              Region 2
            ((C_D_r*(((w3-w2)/(h3-h2))*((4*h3^5-5*h2*h3^4+h2^5)/(20))+(w2*(h3-h2)))))+... Region 3
            ((C_D_c*(((w4-w3)/(h4-h3))*((4*h4^5-5*h3*h4^4+h3^5)/(20))+(w3*(h4-h3)))))+... Region 4
            ((C_D_c*(((w5-w4)/(h5-h4))*((4*h5^5-5*h4*h5^4+h4^5)/(20))+(w4*(h5-h4)))))+... Region 5
            ((C_D_c*(((w6-w5)/(h6-h5))*((4*h6^5-5*h5*h6^4+h5^5)/(20))+(w5*(h6-h5)))))...  Region 6
            )...
            )... 
            / (m*mu*N*(1/L) * cosd(theta2(t)));
            

% Time parameters
t_start = 0.0;
t_end = 1;
dt = 0.001;  % Time step (s)

% Initial conditions
I = 0.0; % Initial value of I(0)
t = t_start;

% Create vectors to store the results for plotting
t_values = t_start:dt:t_end;  % Time values
I_values = zeros(size(t_values));  % Pre-allocate I values

%loop through each value of time and solve for the current at each time
%increment
for i = 1:length(t_values)
    
    %set value of t
    t = t_values(i);
 
    %determine which current function to use depending on t (piecewise):
    if mod(t, 1) < c
        %determine C_D_c (drag coeff. of compliant part of leg - changes
        %depending on direction that the leg is moving in
        if theta1_dot(mod(t,1)) > 0
            C_D_c = 0.5;
        else
            C_D_c = 1;
        end
        %calculate current value for the given t value
        I = I_iso1(mod(t, 1),C_D_c);
    else
        %determine C_D_c (drag coeff. of compliant part of leg - changes
        %depending on direction that the leg is moving in
        if theta2_dot(mod(t,1)) > 0
            C_D_c = 0.5;
        else
            C_D_c = 1;
        end
        %calculate current value for the given t value
        I = I_iso2(mod(t, 1),C_D_c);
    end
        
    %set current I value in the array
    I_values(i) = I;
    
end

% Plot the results
figure;
plot(t_values, I_values, 'LineWidth', 2);
xlabel('t (sec)');
ylabel('I (A)');
title('Solution for Current as a Function of Time');
grid on;

% figure
% hold on
% fplot(theta1)
% fplot(theta2)
% legend("theta1","theta2")

