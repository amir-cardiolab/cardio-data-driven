% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 6 of the paper: 
%Reduced-order physics: Dynamic mode decomposition (DMD) and proper orthogonal decomposition (POD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the code generates the raw data for DMD modes and could be converted to VTK or other formats for visualization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code Courtesy of Mr. Milad Habibi (from Arzani lab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need the optimal_SVHT_coef.m code to run

clc
clear all;
close all;

data1=load ('../data_rsif_paper/VELGround');
Velocity=data1.velocity;
%creating data
% Velocity1=Velocity+normrnd(0,0.2*max(max(Velocity)),size(Velocity));
% Velocity2=Velocity+normrnd(0,0.3*max(max(Velocity)),size(Velocity));
VEL=Velocity;
% Velocity1=zeros(400859,181);
format long
res=zeros(220,1);
res2=zeros(220,1);
qq=size(Velocity);
Velocity1=zeros(qq(1)/3,qq(2));
Velocity2=zeros(qq(1)/3,qq(2));
Velocity3=zeros(qq(1)/3,qq(2));
for i=1:qq(1)
    if rem(i,3)==1
        Velocity1(((i-1)/3)+1,:)=Velocity(i,:);
    end
    if rem(i,3)==2
        Velocity2(((i-2)/3)+1,:)=Velocity(i,:);
    end
    if rem(i,3)==0
        Velocity3(((i)/3),:)=Velocity(i,:);
    end
end
%Velocity=[Velocity1;Velocity2;Velocity3]; 3D calculation
Velocity=[Velocity1;Velocity2];
X1   = Velocity(:,1:end-1);
X2  = Velocity(:,2:end);
[U, S, V] = svd(X1, 'econ');

SIG1=diag(S);
dt=0.95/1500;
beta=(size(S,2)/size(S,1));
thresh1=optimal_SVHT_coef(beta, 0)*median(SIG1);
r =  length(find(diag(S)>thresh1));
r=size(X1,2);
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = U_r' * X2 * V_r / S_r;
[W_r, D] = eig(Atilde);
Phi = X2 * V_r / S_r * W_r; % DMD modes
% aa=vecnorm(Phi);
% Phi=aa.*GramSchmidt(Phi);
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
x1 = X1(:, 1);
b = Phi\x1;
%% DMD reconstruction
mm1 = size(X1, 2)+1; % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1-1)*dt; % time vector
for iter = 1:mm1
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = Phi * time_dynamics;
XDMD2=zeros(size(X1));
XDMD2(:,1)=X1(:,1);
% SPARSE DMD
Fdmd=U_r' * X2 * V_r / S_r;
[Ydmd,Ddmd] = eig(Fdmd);
Edmd = diag(Ddmd);
R = rank(Fdmd)
N=size(V',2);
% Form Vandermonde matrix
Vand = zeros(R,N);
zdmd = Edmd;
for i = 1:N
    
    Vand(:,i) = zdmd.^(i-1);
    
end
L = Ydmd;
R = Vand;
G = S*V';

% Form matrix P, vector q, and scalar s
% J = x'*P*x - q'*x - x'*q + s
% x - optimization variable (i.e., the unknown vector of amplitudes)
P = (L'*L).*conj(R*R');
q = conj(diag(R*G'*L));
s1 = trace(G'*G);

% Cholesky factorization of P
Pl = chol(P,'lower');
% Optimal vector of amplitudes xdmd
amp = (Pl')\(Pl\q);
for i=2:mm1
    XDMD2(:,i)=(X2* V_r / S_r)*(U_r'*XDMD2(:,i-1));
end
ZZ2=((abs(Velocity-real(Xdmd))));
ZZZ2=sum(sum(ZZ2))/(size(Velocity,1)*size(Velocity,2));
PHIB=[];
for i=1:(qq(1)/3)
PHIB=[PHIB;Phi(i,:);Phi((qq(1)/3)+i,:);Velocity3(i,1:end-1)];
end
% cov1=cov(Velocity-Velocity1);
% cov2=cov(Velocity-Velocity2);
% % % cov1=((0.2*max(max(Velocity)))^2)*eye(size(Velocity,1));
% % % cov2=((0.3*max(max(Velocity)))^2)*eye(size(Velocity,1));
% % % [Velocity3]=Kal(XDMD2,Velocity2,Atilde,U_r,cov1,cov2);
% % % [Velocity4,Velocity5]=Kal2(XDMD2,Velocity2,Atilde,U_r,cov2,Velocity1);


% for i=1:size(U_r,2)
%     [c,l]=wavedec(U_r(:,i),1,'db2');
%     si=mad(c)/0.6745;
%     T=si*sqrt(2*log(length(c)));
% end
%     
%data
%XDMD2,Velocity2,U_r,cov1,cov2,