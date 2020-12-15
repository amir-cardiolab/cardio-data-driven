% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 5 of the paper: 
%Data assimilation using the Kalman filter: Merging experimental and computational data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified code from Data-Driven Modeling & Scientific Computation: Methods for Complex Systems & Big Data
%by J. Nathan Kutz (Oxford, 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need Hill_vortex_transient.m 

clear all;

X_IC = [0,0.5,0]; %I.C
dt = 0.001; %delta_t
T_end = 20;
t = 0:dt:T_end ;  %time


%exact solution
[time,X_sol]=ode45('Hill_vortex_transient',t,X_IC);

figure;
plot3(X_sol(:,1),X_sol(:,2),X_sol(:,3));
title('Exact solution');
xlabel('x');
ylabel('y');
zlabel('z');

%noisy solution
Sample_freq = 40; %sampling less frequent also works well.
sample_time = t(1:Sample_freq:end);
n = length(sample_time);
Noise_mu = 0;
Noise_sigma = 0.05;
x_noise = normrnd(Noise_mu,Noise_sigma,[1,n]);
y_noise = normrnd(Noise_mu,Noise_sigma,[1,n]);
z_noise = normrnd(Noise_mu,Noise_sigma,[1,n]);
X_exp = X_sol(1:Sample_freq:end,1) +  x_noise';
Y_exp = X_sol(1:Sample_freq:end,2) + y_noise';
Z_exp = X_sol(1:Sample_freq:end,3) +  z_noise';

figure;
plot3(X_exp,Y_exp,Z_exp);
title('Noisy experimental data');
xlabel('x');
ylabel('y');
zlabel('z');



%Numerical solution with uncertainty
sigma_IC = 0.15  ;
% rng(12345);% for reproducibility of random number generator
% X_IC_p(1) = X_IC(1) + normrnd(0,sigma_IC);  %I.C perturbed
% %rng(1234);% for reproducibility of random number generator
% X_IC_p(2) = X_IC(2) + normrnd(0,sigma_IC);  %I.C perturbed
% %rng(123);% for reproducibility of random number generator
% X_IC_p(3) = X_IC(3) + normrnd(0,sigma_IC);  %I.C perturbed
%perturb with sigma (for sampling consistency)
X_IC_p(1) = X_IC(1) + sigma_IC;
X_IC_p(2) = X_IC(2) + sigma_IC;
X_IC_p(3) = X_IC(3) + sigma_IC;
[time,X_sol_perturbed]=ode45('Hill_vortex_transient',t,X_IC_p);
figure;
plot3(X_sol_perturbed(:,1),X_sol_perturbed(:,2),X_sol_perturbed(:,3));
title('Perturbed solution');
xlabel('x');
ylabel('y');
zlabel('z');


%Kalman filter
x_solution_K = [];
y_solution_K = [];
z_solution_K = [];
for j=1:n-1
    t1 = sample_time(j);
    t2 = sample_time(j+1); 
    t_K = t1:dt:t2;
    [time,X_sol_p]=ode45('Hill_vortex_transient',t_K,X_IC_p);
    xic0=[X_sol_p(end,1); X_sol_p(end,2); X_sol_p(end,3)]; % model estimate
    xdat_exp=[X_exp(j+1); Y_exp(j+1); Z_exp(j+1)]; % experimental data estimate
    K=sigma_IC/(sigma_IC+Noise_sigma); % Kalman gain
    X_IC_p = xic0+(K*(xdat_exp-xic0)); % adjusted state vector
    x_solution_K = [x_solution_K X_sol_p(1:end-1,1)'];
    y_solution_K = [y_solution_K X_sol_p(1:end-1,2)'];
    z_solution_K = [z_solution_K X_sol_p(1:end-1,3)'];
end
 x_solution_K = [x_solution_K X_sol_p(end,1)];
 y_solution_K = [y_solution_K X_sol_p(end,2)];
 z_solution_K = [z_solution_K X_sol_p(end,3)];

figure;
plot3(x_solution_K,y_solution_K,z_solution_K);
title('Kalman solution');
xlabel('x');
ylabel('y');
zlabel('z');


%%%%%%%%%%%% final figures%%%%%

figure;
plot(X_sol(:,1),X_sol(:,2),'color','black','LineWidth',3);
title('Exact solution','FontSize', 24);
xlabel('x','FontSize', 24);
ylabel('y','FontSize', 24);
set(gca,'fontsize',20)

figure;
plot(x_solution_K,y_solution_K,'color','black','LineWidth',3);
title('Kalman solution','FontSize', 24);
xlabel('x','FontSize', 24);
ylabel('y','FontSize', 24);
set(gca,'fontsize',20)

figure;
plot(X_sol_perturbed(:,1),X_sol_perturbed(:,2),'color','black','LineWidth',3);
title('Perturbed solution','FontSize', 24);
xlabel('x','FontSize', 24);
ylabel('y','FontSize', 24);
set(gca,'fontsize',20)




%%%calculate errors %%%%%

error_perturbed = X_sol - X_sol_perturbed;


fprintf('Mean error in Perturbed solution: %f; Max Error: %f\n', mean(sqrt(sum(error_perturbed.^2,2) )), max(sqrt(sum(error_perturbed.^2,2) ) ) );

 X_sol_K = [x_solution_K'  y_solution_K'  z_solution_K'];

error_Kalman = X_sol - X_sol_K;

fprintf('Mean Error in Kalman solution: %f; Max Error: %f\n', mean(sqrt(sum(error_Kalman.^2,2) )), max(sqrt(sum(error_Kalman.^2,2) ))  );

X_exp2 = [X_exp Y_exp Z_exp];  
error_measurement =  X_sol(1:Sample_freq:end,:) - X_exp2 ;

fprintf('Mean Error in measurements (available pts): %f; Max Error: %f\n', mean(sqrt(sum(error_measurement.^2,2) )), max(sqrt(sum(error_measurement.^2,2) ))  );


figure;
semilogy(t,sqrt(sum(error_perturbed.^2,2)),'LineWidth',3);
hold on;
semilogy(t,sqrt( sum(error_Kalman.^2,2) ),'r','LineWidth',3);
xlabel('Time (s)','FontSize', 24);
ylabel('Error','FontSize', 24);
set(gca,'fontsize',20)

