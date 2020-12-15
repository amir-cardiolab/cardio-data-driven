% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 2 of the paper: 
%Principal component analysis (PCA): detecting redundancy and low-dimensionality in data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified code from Data-Driven Science and Engineering: Machine Learning,
%Dynamical Systems, and Control, by Steven L. Brunton and J. Nathan Kutz (Cambridge, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
Flag_data_choice = 0; %if 1 uses Brain aneurysm data. If 0 uses AAA data

if(Flag_data_choice)
    load '../data_rsif_paper/Vel_Ane.mat'; %Brain aneurysm data
else
    load '../data_rsif_paper/AAA_P95.mat';  %Abdominal aortic aneurysm data
end


if(Flag_data_choice)
X = velocity(:,1:4:200);  %downsample to same temporal resolution as AAA (for Vel_Anne)
else
X = velocity(:,1:end);
end

n = size(X,1); %X_nm
m = size(X,2); %X_nm

X_mean = mean(X,2);
%X_std = std(X,0,2);

%subtract mean
for i=1:n  
%X(:,i) = ( X(:,i) - X_mean(i) ) / X_std(i);
X(i,:) = ( X(i,:) - X_mean(i) );
end




X = (1/sqrt(m)) * X;

[u,s,v]=svd(X,'econ'); % perform singular value decomposition (SVD)

sigma = diag(s); %singular values
sigma_energy = cumsum(sigma);



figure;
plot(sigma,'ko','Linewidth',[1.5]);
title('PCA','FontSize', 45);
xlabel('modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

figure;
semilogy(sigma,'ko','Linewidth',[1.5]);
title('PCA, semi-log plot','FontSize', 45);
xlabel('modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
%ylim([1e-1,1e4]);%AA
ylim([1e-2,1e4]);%IA
xlim([0,50]);

figure;
plot(sigma_energy/sum(sigma),'ko','Linewidth',[1.5]);
title('Normalized cummulative energy','FontSize', 45);
xlabel('modes','FontSize', 35);
ylabel('Cummulative energy','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);
