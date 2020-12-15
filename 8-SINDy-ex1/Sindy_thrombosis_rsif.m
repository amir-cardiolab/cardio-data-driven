% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 9 of the paper: 
%Sparse identification of nonlinear dynamics (SINDy): Discovering analytical dynamics
%Example 1: Discover a blood coagulation and thrombosis model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified code from Data-Driven Science and Engineering: Machine Learning,
%Dynamical Systems, and Control, by Steven L. Brunton and J. Nathan Kutz (Cambridge, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need Thrombosis_reduced.m; poolData.m poolDataLIST.m sparsifyDynamics.m 

clear all; 


%% Generate Datak
K = [0.0262; 0.525; 7.223e-6; 5.24e-2; 0.002]; % Papadopoulos et al thrombosis parameters (case 5 in paper)
n = 4;
x0=[0; 9.509e-5; 0; 0.004];  % Initial condition (Papadopoulos 2015 thesis)
tspan=0.:1:5000;
%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n)); %original value
Flag_solver = 0; %if 1 full system, if 0 the reduced thrombosis solver

if(Flag_solver)
options = odeset('RelTol',1e-16,'AbsTol',1e-16*ones(1,n));
[t,x]=ode45(@(t,x) Thrombosis(t,x,K),tspan,x0,options); %full system
else
n=3;
x0_reduced = x0(1:3);
options = odeset('RelTol',1e-16,'AbsTol',1e-16*ones(1,n));
[t,x]=ode45(@(t,x) Thrombosis_reduced(t,x,K,x0),tspan,x0_reduced,options); %reduced to 3 variables
end




figure;
plot(t,x(:,1));
figure;
plot(t,x(:,2));
figure;
plot(t,x(:,3));

figure;
plot3(x(:,1),x(:,2),x(:,3),'Linewidth',3)
xlabel('[IIa]','FontSize', 24);
ylabel('[II]','FontSize', 24);
zlabel('[AP]','FontSize', 24);
set(gca,'fontsize',20)


%% Compute Derivative
%x = x(1:10:end,:); %downsample?

if(Flag_solver)
for i=1:length(x)
    dx(i,:) = Thrombosis(0,x(i,:),K);
end
else
 for i=1:length(x)
    dx(i,:) = Thrombosis_reduced(0,x(i,:),K,x0);
end   
end


%% Build library and compute sparse regression



polyorder = 2; %3
Theta = poolData(x,n,polyorder);  % up to third order polynomials
fprintf('condition number of Theta: %e', cond(Theta));

if(1) %original method
 lambda = 1e-6; %0.0001;  %0.025 this lambda makes too sparse and wrong results;      % lambda is the sparsification knob.
 Xi = sparsifyDynamics(Theta,dx,lambda,n);
end

if(0) %Can compare performance with other optimization methods if desired
   n_d = size(Theta,2);
   cvx_begin;
    variable  Xi(n_d,n);
    minimize( norm( Xi,1) );
    subject to
     Theta* Xi == dx; %this was the original
     %norm(Theta* Xi-dx,1) < 0.01 ; %maybe more relaxed?
   cvx_end;
       
end

if(0)%Can compare performance with other optimization methods if desired
n_d = size(Theta,2);
lambda = 0.01 ;
lambda2 = 0.;
a0 = zeros(n_d,n);
%a0 = Theta\dx;  % initial guess: Least-squares
V = @(a) norm(Theta*a-dx,2) + lambda* norm(a,1) + lambda2 * (norm(a,2)^2);
%opts = psoptimset('MaxIter',19500000);
Xi = patternsearch(V,a0); %global optimization    
end


if(Flag_solver)
poolDataLIST({'[IIa]','[II]','[AP]','[RP]'},Xi,n,polyorder);
else
poolDataLIST({'[IIa]','[II]','[AP]'},Xi,n,polyorder);    
end
