% This Matlab code solves a linear stochastic differential equation (sde) using 
% Euler-Maruyama method. Result is compared with an analytical result, Ref. [1]. 
%
% SDE: dX(t) = -X(t)*dt + dW(t), X(0) = X0, 0 < t < T
%
% It's solution satisfies a following distribution function: 
%                        p(x,t->infty) = (1/sqrt(pi)) * exp(-x^2);
%
% Ref. [1] D. J. Higham, "An algorithm introduction to numerical simulation of stochastic differential equations", 
% SIAM Rev, v43, p525, (2001);
%      
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% Jan 12, 2025 & University of North Dakota
%
function [] = linear_stochastic_diff_eq_2
%
clear; clc;
%
%randn('state',100)
%
T_f = 10000.;   % time t
Nt = 2^16;
dt = T_f/Nt;
%
M = 5.;         % number of paths
%
x0 = rand(M,1);
%
fileID_save_data_1 = fopen('linear_stochastic_diff_eq_2.txt','w');
%
for s = 1:M
    %
    dw = sqrt(dt)*randn(1,Nt);
    for ii = 1:Nt
        %
        dW = sum( dw((ii-1)+1:ii) );
        x = x0 - dt * x0 + dW;
        x0 = x;
        %
        output = [ii, ii*dt, mean(x)];                
        %
        fprintf(fileID_save_data_1, '%4.4f \t %8.4f \t %8.12f\n', output);
        %    
    end
end
%
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_output_data = fopen('linear_stochastic_diff_eq_2.txt', 'r');               
read_output_data = textscan(read_output_data, '%f %f %f');
%t = read_output_data{2};
x_val = read_output_data{3};
%%%


%%% an exact probability distribution function 
x_vals = -5:0.1:5.; 
pdf_exact = (1/sqrt(pi)).*exp(-x_vals.^2);
%%%
nbins = 30;
figure(1)
hold on
plot(x_vals, pdf_exact, 'b-', LineWidth=1.5)   % exact distribution 
histogram(x_val, nbins,'Normalization','pdf')  %  
xlabel('$x$','Interpreter','latex') % ,'fontsize',16
ylabel('Probability distribution function') % , 'Rotation',0 ,'Rotation',1
hold off
set(gca,'FontSize',18)
box on



%%%
return
end
