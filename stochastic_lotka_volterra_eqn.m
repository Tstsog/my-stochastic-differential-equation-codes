% This Matlab code solves the stochastic Lotka-Volterra (stoch LV) model using 
% Euler-Maruyama method [1]. Result is compared with an solution of the deterministic 
% (without stochastic term) equation [2]. 
%
% Simulation of the stochastic Lotka-Volterra model
% dx = (alpha*x - beta*x*y) * dt + sigma_1 * dw1    with x(0) and y(0)
% dy = (-gamma*y + beta*x*y) * dt + sigma_2 * dw2
%
% Ref. [1] D. J. Higham, "An algorithm introduction to numerical simulation of stochastic differential equations", 
% SIAM Rev, v43, p525, (2001);
% Ref. [2] https://en.wikipedia.org/wiki/Lotka–Volterra_equations
%      
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% Jan 23, 2025 & University of North Dakota
%
function [] = stochastic_lotka_volterra_eqn
clear; clc; 
%
alpha = 1.5;
beta = 1.;
delta = 1.;
gamma = 3;
%
sigma_x = 0.1;
sigma_y = 0.3;
%
x0 = 1.5;
y0 = 1.5;
%
ti = 0.;
tf = 25.;
Nt = 5000.;
%
dt = (tf - ti)/Nt;
%
lv_euler(x0, y0, alpha, beta, delta, gamma, Nt, dt);
lv_euler_maruyama(x0, y0, alpha, beta, delta, gamma, sigma_x, sigma_y, Nt, dt);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_output_data = fopen('lv_euler.txt', 'r');               
read_output_data = textscan(read_output_data, '%f %f %f %f');
t_euler = read_output_data{2};
x_euler = read_output_data{3};
y_euler = read_output_data{4};
%
read_output_data = fopen('lv_euler_maruyama.txt', 'r');               
read_output_data = textscan(read_output_data, '%f %f %f %f');
t_euler_maruyama = read_output_data{2};
x_euler_maruyama = read_output_data{3};
y_euler_maruyama = read_output_data{4};

%
figure(1)
hold on
plot(t_euler, x_euler, 'b', 'LineWidth',1.5)
plot(t_euler, y_euler, 'g', 'LineWidth',1.5)
%
plot(t_euler_maruyama, x_euler_maruyama, 'b--', 'LineWidth',1.5)
plot(t_euler_maruyama, y_euler_maruyama, 'g--', 'LineWidth',1.5)
hold off
xlabel('$time$', 'interpreter','latex')
ylabel('$population$','interpreter','latex')
set(gca,'FontSize',18)
box on

%%%
figure(2)
hold on
plot(x_euler, y_euler, 'b', 'LineWidth',1.5)
plot(x_euler_maruyama, y_euler_maruyama, 'b--', 'LineWidth',1.5)
hold off
xlabel('$x(t)$','interpreter','latex')
ylabel('$y(t)$','interpreter','latex')
set(gca,'FontSize',18)
box on


%%%
return
end

%%%
function lv_euler(x0, y0, alpha, beta, delta, gamma, Nt, dt)
%
fileID_save_data_1 = fopen('lv_euler.txt','w');
%
for ii = 1:Nt
    %
    x = x0 + dt * (alpha * x0 - beta * x0 * y0);
    y = y0 + dt * (delta * x0 * y0 - gamma * y0);
    %
    x0 = x;
    y0 = y;
    %
    output = [ii, ii*dt, x, y];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.4f \t %8.12f \t %8.12f\n', output);
    %     
end
%
fclose(fileID_save_data_1);
%%%
return
end

%%%
function lv_euler_maruyama(x0, y0, alpha, beta, delta, gamma, sigma_x, sigma_y, Nt, dt)
%
dw1 = sqrt(dt)*randn(1,Nt);                      % Brownian increments
dw2 = sqrt(dt)*randn(1,Nt);                      % Brownian increments
%
fileID_save_data_1 = fopen('lv_euler_maruyama.txt','w');
%
for ii = 1:Nt
    %
    dW1_inc = sum(dw1((ii-1)+1:ii));                % dW1(t) = W1(i) - W1(i-1)
    dW2_inc = sum(dw2((ii-1)+1:ii));                % dW2(t) = W2(i) - W2(i-1)        
    %
    x = x0 + dt * (alpha * x0 - beta * x0 * y0) + sigma_x * x0 * dW1_inc;
    y = y0 + dt * (delta * x0 * y0 - gamma * y0) + sigma_y * y0 * dW2_inc;
    %
    x0 = x;
    y0 = y;
    %
    output = [ii, ii*dt, x, y];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.4f \t %8.12f \t %8.12f\n', output);
    %     
end
%
fclose(fileID_save_data_1);
%%%
return
end
