function myRPCA_rsif
% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 3 of the paper: 
%Robust principal component analysis (RPCA): noisy and fluctuating data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified code from Data-Driven Science and Engineering: Machine Learning,
%Dynamical Systems, and Control, by Steven L. Brunton and J. Nathan Kutz (Cambridge, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



load '../data_rsif_paper/AAA_P95.mat'  %AAA velocity data


X = velocity(:,1:end);

X_z = X(3:3:end,:);
Noise_level = 0.1*max(abs(X_z(:))) ; % SD of the Gaussian noise added
Fraction_data_noise = 0.3; %fraction of data where noise is added





%Add noise 
[Xn,Yn]=size(X);
n_pts = Xn*Yn;


n_noise = floor(Fraction_data_noise*n_pts);    %number of data pts for adding noise 
r1 = randintrlv(1:n_pts,1234); %random noise location
samples_noise = r1(1:n_noise);
X_temp = reshape(X,n_pts,1);
X_temp(samples_noise) = X_temp(samples_noise) + Noise_level*randn(size(X_temp(samples_noise)));   
X_noise = reshape(X_temp,Xn,Yn);



X_mean = mean(X,2);
%X_std = std(X,2);
%subtract mean 
X2 = X;
for i=1:Xn  
X2(i,:) = ( X2(i,:) - X_mean(i) ) ;
end
[u,s,v]=svd((1/sqrt(Yn))*X2,'econ'); % perform singular value decomposition (SVD)


X_mean = mean(X_noise,2);
%subtract mean 
X_noise2 = X_noise;
for i=1:Xn  
X_noise2(i,:) = ( X_noise2(i,:) - X_mean(i) );
end
[u_n,s_n,v_n]=svd( (1/sqrt(Yn))*X_noise2,'econ');


%%%% RPCA on original data%%%%%%
mu = Xn*Yn/(4*sum(abs(X(:))));
lambda = 1/sqrt(max(Xn,Yn));
thresh = 1e-6*norm(X,'fro');
L = zeros(size(X));
S = zeros(size(X));
Y = zeros(size(X));
count = 0;
while((norm(X-L-S,'fro')>thresh) && (count<1500))  % count<1000
L = SVT(X-S+(1/mu)*Y,1/mu);
S = shrink(X-L+(1/mu)*Y,lambda/mu);
Y = Y + mu*(X-L-S);
count = count + 1;
count
norm(X-L-S,'fro')
end

L_mean = mean(L,2);
%subtract mean
for i=1:Xn  
L(i,:) = ( L(i,:) - L_mean(i) ) ;
end
[u_L,s_L,v_L]=svd((1/sqrt(Yn))*L,'econ');



%%%% RPCA on noisy data%%%%%%
mu = Xn*Yn/(4*sum(abs(X_noise(:))));
lambda = 1/sqrt(max(Xn,Yn));
thresh = 1e-6*norm(X,'fro');
L = zeros(size(X_noise));
S = zeros(size(X_noise));
Y = zeros(size(X_noise));
count = 0;
while((norm(X_noise-L-S,'fro')>thresh) && (count<1500))  % count<1000
L = SVT(X_noise-S+(1/mu)*Y,1/mu);
S = shrink(X_noise-L+(1/mu)*Y,lambda/mu);
Y = Y + mu*(X_noise-L-S);
count = count + 1;
count
norm(X_noise-L-S,'fro')
end

L_mean = mean(L,2);
%subtract mean 
for i=1:Xn  
L(i,:) = ( L(i,:) - L_mean(i) ) ;
end
[u_Ln,s_Ln,v_Ln]=svd((1/sqrt(Yn))*L,'econ');


%%%% plot%%%%%
sigma = diag(s); %singular values, original data
sigma_energy = cumsum(sigma);

sigma_n = diag(s_n); %singular values, noisy data
sigma_energy_n = cumsum(sigma_n);

sigma_L = diag(s_L); %singular values, RPCA on original data
sigma_energy_L = cumsum(sigma_L);

sigma_Ln = diag(s_Ln); %singular values, RPCA on noisy data
sigma_energy_Ln = cumsum(sigma_Ln);


%%% Original data
figure;
plot(sigma,'ko','Linewidth',[1.5]);
title('Principal components','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

figure;
semilogy(sigma,'ko','Linewidth',[1.5]);
title('Principal components, semi-log plot','FontSize', 34);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);
ylim([1e0,1e4]);

figure;
plot(sigma_energy/sum(sigma),'ko','Linewidth',[1.5]);
title('Normalized cummulative energy','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Cummulative energy','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

%%%%% Noisy data
figure;
plot(sigma_n,'ko','Linewidth',[1.5]);
title('Principal components','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

figure;
semilogy(sigma_n,'ko','Linewidth',[1.5]);
title('Principal components, semi-log plot','FontSize',34);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);
ylim([1e0,1e4]);

figure;
plot(sigma_energy_n/sum(sigma_n),'ko','Linewidth',[1.5]);
title('Normalized cummulative energy','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Cummulative energy','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

%%%%
%%%%% RPCA on original data
figure;
plot(sigma_L,'ko','Linewidth',[1.5]);
title('Principal components','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

figure;
semilogy(sigma_L,'ko','Linewidth',[1.5]);
title('Principal components, semi-log plot','FontSize', 34);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);
ylim([1e-14,1e4]);

figure;
plot(sigma_energy_L/sum(sigma_L),'ko','Linewidth',[1.5]);
title('Normalized cummulative energy','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Cummulative energy','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

%%%%% RPCA on noisy data
figure;
plot(sigma_Ln,'ko','Linewidth',[1.5]);
title('Principal components','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

figure;
semilogy(sigma_Ln,'ko','Linewidth',[1.5]);
title('Principal components, semi-log plot','FontSize', 34);
xlabel('Modes','FontSize', 35);
ylabel('Singular values','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);
ylim([1e-14,1e4]);


figure;
plot(sigma_energy_Ln/sum(sigma_Ln),'ko','Linewidth',[1.5]);
title('Normalized cummulative energy','FontSize', 45);
xlabel('Modes','FontSize', 35);
ylabel('Cummulative energy','FontSize', 35);
set(gca,'fontsize',29)
xlim([0,50]);

%save all data
%save('myRPCAresults.mat')

end

function out = shrink(X,tau)
 out = sign(X).*max(abs(X)-tau,0);
end

function out = SVT(X,tau)
[U,S,V] = svd(X,'econ');
out = U*shrink(S,tau)*V';
end






