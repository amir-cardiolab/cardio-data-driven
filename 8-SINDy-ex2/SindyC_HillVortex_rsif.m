% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 9 of the paper: 
%Sparse identification of nonlinear dynamics (SINDy): Discovering analytical dynamics
%Example 2: Discover analytical velocity in transient vortical flows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified code from Data-Driven Science and Engineering: Machine Learning,
%Dynamical Systems, and Control, by Steven L. Brunton and J. Nathan Kutz (Cambridge, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need Hill_vortex_transient.m; poolData_control.m poolDataLIST_unsteady.m sparsifyDynamics.m 



clear all; close all; 

%% Generate Data

%X_IC = [0.02,0.5,.01]; %I.C
X_IC = [0.02,0.05,.01]; %I.C
dt = 0.001; %delta_t
T_end = 8;
t = 0:dt:T_end ;  %time
n = 3;
[t,x]=ode45('Hill_vortex_transient',t,X_IC);


figure;
plot(t,x(:,1));
figure;
plot(t,x(:,2));
figure;
plot(t,x(:,3));

figure;
plot3(x(:,1),x(:,2),x(:,3),'r','Linewidth',3)
xlabel('x','FontSize', 24);
ylabel('y','FontSize', 24);
zlabel('z','FontSize', 24);
set(gca,'fontsize',20)

%% Compute Derivative
%x = x(1:10:end,:); %downsample?
for i=1:length(x)
    dx(i,:) = Hill_vortex_transient(t(i),x(i,:));
end


%% Build library and compute sparse regression
polyorder = 2; %3
Theta = poolData_control(x,n,polyorder,t);  % up to third order polynomials
if(1) %original method
 lambda = 0.01; % lambda is the sparsification knob.
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
lambda = 0.001 ;
lambda2 = 0.;
a0 = zeros(n_d,n);
%a0 = Theta\dx;  % initial guess: Least-squares
V = @(a) norm(Theta*a-dx,2) + lambda* norm(a,1) + lambda2 * (norm(a,2)^2);
%opts = psoptimset('MaxIter',19500000);
Xi = patternsearch(V,a0); %global optimization    
end


poolDataLIST_unsteady({'x','y','z'},Xi,n,polyorder);

