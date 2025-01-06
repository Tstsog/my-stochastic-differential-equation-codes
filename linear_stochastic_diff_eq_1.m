% This Matlab code solves a linear stochastic differential equation (sde) using 
% Euler-Maruyama method. Result is compared with an analytical result, Ref. [1]. 
%
% SDE: dX(t) = a*X(t)*dt + b*X(t)*dW(t), X(0) = X0, 0 < t < T
%
% analytic result: X(t) = X(0) * exp((a - 0.5*b^2)*t + b*W(t))
%
% Ref. [1] D. J. Higham, "An algorithm introduction to numerical simulation of stochastic differential equations", 
% SIAM Rev, v43, p525, (2001);
%      
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% Jan 6, 2025 & University of North Dakota
%
function [] = linear_stochastic_diff_eq_1
%
clear; clc;
%
randn('state',100)
%
a = 0.25;
b = 1.20;
x0 = 2.;
%
T = 2.5;
Nt = 2^10;
dt = T/Nt;
%
dw = sqrt(dt)*randn(1,Nt);                      % Brownian increments 
w = cumsum(dw);
%
time = dt:dt:T;
x_exact = x0 * exp((a - 0.5*b^2)*time + b * w);  % exact analytical formula
%
fileID_save_data_1 = fopen('linear_stochastic_diff_eq_1.txt','w');
%
for ii = 1:Nt
    %
    dW = sum(dw((ii-1)+1:ii));                % dW(t) = W(i) - W(i-1)
    x = x0 + dt * a * x0 + b * x0 * dW;
    x0 = x;
    %
    output = [ii, ii*dt, x];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.4f \t %8.12f\n', output);
    %    
end
%
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_output_data = fopen('linear_stochastic_diff_eq_1.txt', 'r');               
read_output_data = textscan(read_output_data, '%f %f %f');
t = read_output_data{2};
x_val = read_output_data{3};
%%%
figure(1)
hold on
plot([0,time(1:5:length(time))], [x0,x_exact(1:5:length(x_exact))], 'ro', 'LineWidth',1.5)
plot([0;t], [x0; x_val], 'b', 'LineWidth',1.5 )
hold off
xlabel('\mbox{Time}','Interpreter','latex') % ,'fontsize',16
ylabel('$X(t)$','Interpreter','latex', 'Rotation',1) % , 'Rotation',0
%axis([0. 8. -2.860 -2.845])
set(gca,'FontSize',18)
box on


%%%
return
end
