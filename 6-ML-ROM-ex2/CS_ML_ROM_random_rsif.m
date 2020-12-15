% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 7 of the paper: 
%Machine learning reduced-order models (ROM): Overcoming uncertainty in
%computational models using low-fidelity experimental data
%Example 2: Reconstruct high-resolution cerebral aneurysm flow with uncertain viscosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need the cosamp.m function to run
%Also, To run with the cvx_begin method need to install:
%CVX: a Matlab-based convex modeling framework
%http://cvxr.com/cvx/
%%%%%%%%%%%%%%%%%
% Compressed sensing + ML 
%This code is based on random selection of data + noise added
%CVX method seems to work much better than cosamp here (it can identify
%which modes are active so it can be used for parameter identification)
%%%%%%%%%%%%%%%%
%final parameters tested:
%nSensors = 300 (This is really 2* number of sensors; 2D vector data)
%NoiseLevel = 0.01;,eps = 0.1+10 and L1 norm worked.
%NoiseLevel=0, eps=0.1 and L2 norm worked
%Noise Level=0.1, eps = 0.1+90; and L1 works (changing eps to 80 seems to
%cause solution to be not sparse); these give poor reconstruction error,
%but manage to identify the parameter; if we slowly decrease the number of
%sensors we see that we get poor reconstruction but we still can identify
%from the most dominant mode (other modes start to become more dominant)
%L2 norm, eps=10 and Noise=0.05 (L1 does not work here!)


clear all;

data_file = '../data_rsif_paper/IA_mu_data/vel_mu';
mu = linspace(0.03,0.2,8); %range of viscosity values

p=8; %number of parameters
n_modes = 12; %20; % number of modes to keep (could optimize this for every parameter)






%load data and do POD

k=1;
for i=0:p-1
    file_name2 = sprintf('%s%i.mat',data_file,i);   
    data =load (file_name2,'Vel_2D_unsteady');
    data = data.Vel_2D_unsteady;
    if (i==0)
        n_data_pts = size(data,1);
        Psi = zeros(n_data_pts,p*n_modes);
        %data_all = zeros(n_data_pts, size(data,2)*p );
        data_all=[];
        %DataMean_all = [];
        data_all_meansub = [];
        Psi = [];
    end
    data_all = [data_all, data];
    % subtract mean
    DataMean = mean(data,2);
    %DataMean_all = [DataMean_all, DataMean];
    nSnaps = size(data,2);
    DataSub = zeros( n_data_pts , nSnaps );
    for j=1:nSnaps
     DataSub(:,j) = data(:,j); % - DataMean; 
    end
    data_all_meansub = [data_all_meansub, DataSub];
    [U,S,V] = svd(DataSub,'econ');
    %Psi(:,k:k+n_modes-1) = U(:,1:n_modes);
    %k = k+n_modes;
    Psi = [Psi, U(:,1:n_modes)];
end




snapind = 390 ; %60 % snapshot to reconstruct
nSensors = 300; %it is 2D vel data so 300/2 sensors %Should be less than the total radial points (herein, 7886)
nModesSparse = p * n_modes ; % cheating a bit by using a reduced basis here
NoiseLevel = 0.05;
rng(12345); % for reproducable results
%perm = round(rand(nSensors, 1) * n); % choose random sensor locations
r1 = randintrlv(1:n_data_pts,12345);
perm = r1(1:nSensors); % choose random sensor locations
%y = data_all(perm,snapind); % compressed measurement
y = data_all_meansub(perm,snapind); % compressed measurement
%y = y + NoiseLevel*randn(size(y)); % add noise, if desired
y = y + NoiseLevel*randn(size(y)).*y; % add noise (relative to data)
Theta = Psi(perm,1:nModesSparse);
s = cosamp(Theta,y,floor(nModesSparse/2),1.e-10,100); %reconstruct POD coefficients using L1 signal recovery method
%s = cosamp(Theta,y,n_modes,1.e-10,1000);

Error1 = abs( Psi(:,1:nModesSparse)*s - (data_all(:,snapind) ));
fprintf('Error in random sensors = %f\n', norm(Error1) );
fprintf('Relative error in random sensors = %f\n', norm(Error1) / mean ( abs(data_all(:,snapind) )) );


%%cvx optimization (the constraint needs to be exactly satisfied!!)
%p=8;
eps = 7; 
cvx_begin;
variable s(nModesSparse);
 minimize( norm(s,1) );
 subject to
  %Theta*s == y;
  norm(Theta*s - y,2) <= eps ; 
cvx_end;

Error1 = abs( Psi(:,1:nModesSparse)*s - (data_all(:,snapind) ));
fprintf('cvx: Error in random sensors = %f\n', norm(Error1) );
fprintf('cvx: Relative error in random sensors = %f\n', norm(Error1) / mean ( abs(data_all(:,snapind) )) );

%Save the results to text
soln = Psi(:,1:nModesSparse)*s;
soln = reshape(soln,2,n_data_pts/2);
soln = soln';
fid = fopen('random-sensors.txt','wt');
fprintf(fid,'%d %d\n',soln.');
fclose(fid);
%save the original data to text
soln = data_all(:,snapind);
soln = reshape(soln,2,n_data_pts/2);
soln = soln';
fid = fopen('original-data.txt','wt');
fprintf(fid,'%d %d\n',soln.');
fclose(fid);
 

figure;
x = 1:p * n_modes;
plot(x ,s,'*','Linewidth',2);
hold on;
xlabel('Modes','FontSize', 24);
ylabel('Mode coefficient','FontSize', 24);
title('Mode identification (random sensors)','FontSize', 24);
set(gca,'fontsize',20)
set(gcf,'Position',[100 100 1400 300]);
%C = {'k','b','r','g','y','m','c',[1 0.4 0.6] }; % Cell array of colrors (for p=8)

Yl=ylim;
y1 = [Yl(1) Yl(1) Yl(2) Yl(2) ];

xv = n_modes;  
txt_pos = xv/2-4;
for i=1:p
 x1=[0 xv xv 0];   
 patch(x1,y1,'r','FaceAlpha',0.05);
 mymu=sprintf('%s=%.2f','\mu',mu(i));
 text( txt_pos,(Yl(1) + Yl(2))/2,mymu,'FontSize',20)
 xv = xv + n_modes;
 txt_pos = txt_pos + n_modes;
end



