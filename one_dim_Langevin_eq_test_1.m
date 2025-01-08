% This Matlab code solves a one-dimensional (one-dime) Langevin equation (eq) using 
% Euler-Maruyama method Ref. [1]. Average of velocity squared (v^2) is obtained as well.  
%
% Langevin equation : % dx/dt = v(t)
%                       m*dv(t)/dt = -gamma*v(t) + xi(t), X(0) = X0, 0 < t < T
%
% Ref. [1] D. J. Higham, "An algorithm introduction to numerical simulation of stochastic differential equations", 
% SIAM Rev, v43, p525, (2001);
%      
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% Jan 6, 2025 & University of North Dakota
%
function [] = one_dim_Langevin_eq_test_1
%
clear; clc;
%
%randn('state',100)
%
m = 2.5;
gamma = 1.5;
%
T_f = 10000.;
Nt = 2^18;
dt = T_f/Nt;
%
x0 = randn(1);
v0 = 0.;
%
dw = sqrt(dt)*randn(1,Nt);                   % Brownian increments 
%
fileID_save_data_1 = fopen('one_dim_Langevin_eq_test_1.txt','w');
%
for ii = 1:Nt
    %
    dW = sum(dw((ii-1)+1:ii));                % dW(t) = W(i) - W(i-1)
    x = x0 + dt * v0;
    v = v0 + dt * (-gamma/m)*v0 + (1/m) * dW;
    %
    x0 = x;
    v0 = v;
    %
    output = [ii, ii*dt, x, v.*v];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.4f \t %8.12f \t %8.12f\n', output);
    %    
end
%
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_output_data = fopen('one_dim_Langevin_eq_test_1.txt', 'r');               
read_output_data = textscan(read_output_data, '%f %f %f %f');
t = read_output_data{2};
x_val = read_output_data{3};
v_val = read_output_data{4};
%%%

%%%
% sigma = 1/2
average_of_v_squared = sum(v_val)/Nt % <0.5*m*v^2> = sigma/(2*gamma)
%average_of_v_squared = 0.1328 vs exact value = 0.1333...

%
figure(1)
hold on
plot(t, x_val, 'b', 'LineWidth',1.5 )
%plot(t, v_val, 'g', 'LineWidth',1.5 )
hold off
xlabel('\mbox{Time}','Interpreter','latex') % ,'fontsize',16
ylabel('$x(t)$','Interpreter','latex') % , 'Rotation',0
%axis([0. 8. -2.860 -2.845])
set(gca,'FontSize',18)
box on

%%%
return
end
