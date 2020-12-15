% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 6 of the paper: 
%Reduced-order physics: Dynamic mode decomposition (DMD) and proper orthogonal decomposition (POD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the code generates the raw data for POD modes and could be converted to VTK or other formats for visualization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code Courtesy of Mr. Milad Habibi (from Arzani lab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all;
close all;

data1=load ('../data_rsif_paper/VELGround');
Velocity=data1.velocity;
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
X=Velocity;
C=X'*X;
[Q1,Omega1]=eig(C);
[d,ind] = sort(diag(Omega1),'descend');
Omega2=Omega1(ind,ind);
Omega=(Omega1(ind,ind)^-0.5);
% Omega=Omega(1:150,:);
Q=Q1(:,ind);
% Q=Q(:,1:150);
Phi=X*Q*(Omega);
B=Phi'*X;
Xt=Phi*B;
ZZ2=((abs(Xt-real(X))));
ZZZ2=sum(sum(ZZ2))/(qq(1)*qq(2));
PHIB=[];
for i=1:(qq(1)/3)
PHIB=[PHIB;Phi(i,:);Phi((qq(1)/3)+i,:);Velocity3(i,:)];
end




