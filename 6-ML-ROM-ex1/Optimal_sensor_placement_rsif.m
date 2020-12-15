% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 7 of the paper: 
%Machine learning reduced-order models (ROM): Overcoming uncertainty in
%computational models using low-fidelity experimental data
%Example 1: Reconstruct high-resolution Womersley flow with uncertain Womersley number and optimal sensor placement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need the cosamp.m function to run

clear all
close all
clc

load ../data_rsif_paper/Vel_Wom_alpha_ranges_complete.mat% Wo=1-29 %V_data_all and radius arrays
snapind = 17*100 + 1; % snapshot to reconstruct => this is Wo=14.68 ; 
n = size(V_data_all,1);
nSnaps = size(V_data_all,2);
%% Do SVD

% compute mean
DataMean = mean(V_data_all,2);
% subtract mean
for j=1:nSnaps
DataSub(:,j) = V_data_all(:,j) - DataMean; 
end



[U,S,V] = svd(DataSub,'econ');
%% Plot left singular vectors (POD modes)
figure;
for ii = 1:8
    subplot(4,2,ii)
    plot(radius, U(:,ii));
    title(['Mode ',num2str(ii)])
end



%% Try compressed sensing using random sensor locations

if(0) %plot the energy content to decide on the number of modes
    sigma = diag(S); %singular values
    %sigma_energy = cumsum(sigma);
    figure;
    plot(sigma,'ko');
    title('principal components');
end

    


nSensors = 5; %Should be less than the total radial points (herein, 150)
nModesSparse = 10; % cheating a bit by using a reduced basis here
NoiseLevel = 0.;
Random_number = 1367;  % for reproducable results
rng(Random_number); % for reproducable results
%perm = round(rand(nSensors, 1) * n); % choose random sensor locations
r1 = randintrlv(1:n,Random_number);
perm = r1(1:nSensors); % choose random sensor locations
y = DataSub(perm,snapind); % compressed measurement
y = y + NoiseLevel*randn(size(y)); % add noise, if desired
Theta = U(perm,1:nModesSparse);
s = cosamp(Theta,y,floor(nModesSparse/2),1.e-10,100); %reconstruct POD coefficients using L1 signal recovery method

%%





figure;
plot(radius,DataMean+DataSub(:,snapind),'Linewidth',3);
hold on;
plot(radius(perm),DataMean(perm)+DataSub(perm,snapind),'ko','Linewidth',10)
title('True data and random sensor locations','FontSize', 25)   
xlabel('Radius (m)','FontSize', 25);
ylabel('Velocity (m/s)','FontSize', 25);
set(gca,'fontsize',20)
xlim([0,.012]);

figure;
plot(radius,DataMean+U(:,1:nModesSparse)*s,'Linewidth',3 );
title('Reconstructed; random sensors with clean data','FontSize', 25)   
xlabel('Radius (m)','FontSize', 25);
ylabel('Velocity (m/s)','FontSize', 25);
set(gca,'fontsize',20)
xlim([0,.012]);



Error1 = abs( U(:,1:nModesSparse)*s - DataSub(:,snapind) );
fprintf('Error in random sensors with clean data= %f\n', norm(Error1) );


%% Do optimal sensor placement version
r = 5; % 10 % r optimally-placed sensors
[~,~,pivot] = qr(U(:,1:r)','vector');
sensors = pivot(1:r);

% this is a = (C Phi_r)^-1 y:
y = DataSub(sensors,snapind);
y = y + NoiseLevel*randn(size(y));
% reconstruct state in low-dimensional space
a = U(sensors,1:r)\y;

% Estimate of full state
%this is x = Phi_r a
DataReconLS = U(:,1:r)*a;



figure;
plot(radius,DataMean+DataSub(:,snapind),'Linewidth',3);
hold on;
plot(radius(sensors),DataMean(sensors)+DataSub(sensors,snapind),'ko','Linewidth',10)
title('True data and optimal sensor locations','FontSize', 25)
xlabel('Radius (m)','FontSize', 25);
ylabel('Velocity (m/s)','FontSize', 25);
set(gca,'fontsize',20)
xlim([0,.012]);

figure;
plot(radius,DataMean+DataReconLS,'Linewidth',3 );
title('Reconstructed; optimal sensors with clean data','FontSize', 25);  
xlabel('Radius (m)','FontSize', 25);
ylabel('Velocity (m/s)','FontSize', 25);
set(gca,'fontsize',20)
xlim([0,.012]);



Error2 = abs ( DataReconLS - DataSub(:,snapind) ) ;
fprintf('Error in optimal sensors with clean data= %f\n', norm(Error2) );



%% Compare condition numbers
% it is the large condition number of Theta that contributes to the
% reconstruction using random sensors being sensitive to noise

%cond(Theta) % for random placement of sensors
%cond(U(sensors,1:r)) % for optimal placement

%% Now repeat with noisy data


Max_vel = max( abs( DataMean+DataSub(:,snapind) ) );
NoiseLevel_perc = 0.1; %Noise (fraction of max velocity)
NoiseLevel = NoiseLevel_perc * Max_vel;
rng(Random_number); % for reproducable results
%perm = round(rand(nSensors, 1) * n); % choose random sensor locations
r1 = randintrlv(1:n,Random_number);
perm = r1(1:nSensors); % choose random sensor locations
y = DataSub(perm,snapind); % compressed measurement
y = y + NoiseLevel*randn(size(y)); % add noise, if desired
Theta = U(perm,1:nModesSparse);
s = cosamp(Theta,y,floor(nModesSparse/2),1.e-10,100);

%%



figure;
plot(radius,DataMean+U(:,1:nModesSparse)*s,'Linewidth',3 );
title('Reconstructed; random sensors with noisy data','FontSize', 25)
xlabel('Radius (m)','FontSize', 25);
ylabel('Velocity (m/s)','FontSize', 25);
set(gca,'fontsize',20)
xlim([0,.012]);    



Error3 = abs( U(:,1:nModesSparse)*s - DataSub(:,snapind) );
fprintf('Error in random sensors with noisy data= %f\n', norm(Error3) );

%% Do optimal sensor placement version
%r = 10;
[~,~,pivot] = qr(U(:,1:r)','vector');
sensors = pivot(1:r);

% this is a = (C Phi_r)^-1 y:
y = DataSub(sensors,snapind);
y = y + NoiseLevel*randn(size(y));
% reconstruct state in low-dimensional space
a = U(sensors,1:r)\y;

% Estimate of full state
%this is x = Phi_r a
DataReconLS = U(:,1:r)*a;

%%




figure;
plot(radius,DataMean+DataReconLS,'Linewidth',3 );
title('Reconstructed; optimal sensors with noisy data','FontSize', 25)
xlabel('Radius (m)','FontSize', 25);
ylabel('Velocity (m/s)','FontSize', 25);
set(gca,'fontsize',20)
xlim([0,.012]);    





Error4 = abs ( DataReconLS - DataSub(:,snapind) ) ;
fprintf('Error in optimal sensors with noisy data= %f\n', norm(Error4) );







