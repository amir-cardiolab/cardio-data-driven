% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 8 of the paper: 
%Low-rank data recovery from random spatiotemporal measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!Requires Two packages to be installed first:
%Sparco Toolbox 
% http://www.cs.ubc.ca/labs/scl/sparco/
%Also, need to install spot-master
%https://github.com/mpf/spot
%Need need to first run addpath /directory-to/spot-master
%Need IST_MC.m

clear all;
Flag_data_choice = 1; %If 1 uses the 2D aneurysm model; if 0 uses the 3D coronary stenosis model
Frac_sample = .4; %e.g., 0.4 or 0.1 %fraction of the data sampled randomly

if (Flag_data_choice)
load '../data_rsif_paper/vel_2Daneu_crop.mat' %velocity data for the 2D aneurysm model in the region of interest
else
load '../data_rsif_paper/Vel_Ste_51times.mat' %velocity data for the 3D coronary stenosis model   
end





% Vectrozizing the matrix for projection
x = velocity(:);

sizeX = size(velocity);

% Generating random sampling points
T = randperm(prod(sizeX));
IDX = T(1:round(Frac_sample*prod(sizeX))); %sampling
% Creating operator for selecting entries at the chosen random locations 
% Requires Sparco Toolbox 
% http://www.cs.ubc.ca/labs/scl/sparco/
A = opRestriction(prod(sizeX), IDX);

% Sampled data
%y = A(x,1);
y = A*x;




XRec = IST_MC(y,A,sizeX); % Regularized Iterated Soft Thresholding




%%%%%%%%%%%%%%%%%%% Final plots %%%%%%%%%


 
 %compare temporal plot at a point

if (Flag_data_choice) 
t = linspace(0,0.95,size(velocity,2) );
     
 n=2477; %n should be odd such that n+1 is the same point
vel_vec = velocity(n,:);  
vel_vec2 = velocity(n+1,:); 
XRec_vec = XRec(n,:); 
XRec_vec2 = XRec(n+1,:); 

vel_mag = sqrt(  vel_vec.^2 + vel_vec2.^2   );
XRec_mag = sqrt( XRec_vec.^2 + XRec_vec2.^2   );
figure; plot(t,vel_mag,'Linewidth',3); hold on; plot(t,XRec_mag,'r','Linewidth',3)
xlabel('Time (s)','FontSize', 24);
ylabel('Velocity (cm/s)','FontSize', 24);
set(gca,'fontsize',20)

n = 789; %n should be odd such that n+1 is the same point
vel_vec = velocity(n,:);  
vel_vec2 = velocity(n+1,:); 
XRec_vec = XRec(n,:); 
XRec_vec2 = XRec(n+1,:); 

vel_mag = sqrt(  vel_vec.^2 + vel_vec2.^2   );
XRec_mag = sqrt( XRec_vec.^2 + XRec_vec2.^2   );
figure; plot(t,vel_mag,'Linewidth',3); hold on; plot(t,XRec_mag,'r','Linewidth',3)
xlabel('Time (s)','FontSize', 24);
ylabel('Velocity (cm/s)','FontSize', 24);
set(gca,'fontsize',20)


else

t = linspace(0,1,size(velocity,2) );
    
%Reconstructed Point: [-10.103	7.1791	-5.2222]
n = 15022; %3D data: n should be such that n+1, n+2 is the same point (1,4,7...)
vel_vec = velocity(n,:);  
vel_vec2 = velocity(n+1,:); 
XRec_vec = XRec(n,:); 
XRec_vec2 = XRec(n+1,:); 

vel_mag = sqrt(  vel_vec.^2 + vel_vec2.^2   );
XRec_mag = sqrt( XRec_vec.^2 + XRec_vec2.^2   );
figure; plot(t,vel_mag,'Linewidth',3); hold on; plot(t,XRec_mag,'r','Linewidth',3)
xlabel('Time (s)','FontSize', 24);
ylabel('Velocity (cm/s)','FontSize', 24);
set(gca,'fontsize',20)

end
 
 
 