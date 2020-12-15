% Code used in "Data-driven cardiovascular flow modeling: examples and
% opportunities" by Arzani & Dawson.
%Paper: https://arxiv.org/abs/2010.00131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 4 of the paper: 
%Compressed sensing (CS): Reconstructing high-resolution data from low-resolution sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified code from Data-Driven Science and Engineering: Machine Learning,
%Dynamical Systems, and Control, by Steven L. Brunton and J. Nathan Kutz (Cambridge, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!! To run with the cvx_begin method need to install:
%CVX: a Matlab-based convex modeling framework
%http://cvxr.com/cvx/

clear all;

Flag_transform = 1; %if 1 DCT; if 0 FFT
Flag_remove_DC = 1; %if one removes the DC component of signal (better CS performance)



%Patient51 AAA PCMRI data infrarenal aorta waveform.
n=1000;
T=20;
t = linspace(0,T,n); %time 
        a0 =   18.4 ;
       a1 =       15.98  ;
       b1 =       27.19 ; 
       a2 =       -26.7 ; 
       b2 =       18.07 ;
       a3 =      -4.011 ; 
       b3 =      -6.185;  
       w =       6.198;
Vel = a0 + a1*cos(t*w) + b1*sin(t*w) + a2*cos(2*t*w) + b2*sin(2*t*w) + a3*cos(3*t*w) + b3*sin(3*t*w);
    


if (Flag_remove_DC)
    Mean_vel = mean(Vel);
    Vel = Vel - Mean_vel ;
end



m=100;
r1 = randintrlv(1:n,1234);
samples = r1(1:m);


Vel_coarse = Vel(samples);
max_vel = max(Vel_coarse);

% Compressed sensing problem
if (Flag_transform==1)
 Psi = dct(eye(n, n)); 
else
 Psi = fft(eye(n, n));   
end

Theta = Psi(samples, :); % Measure rows of Psi

%We can try different optimization approaches here:
if(0)
 s = cosamp(Theta,Vel_coarse',5,1e-10,100); % compressed sensing via matching pursuit
else if(1) %This is the one used in the paper
   cvx_begin;
    variable s(n);
    minimize( norm(s,1) );
    subject to
     Theta*s == Vel_coarse'; %this was the original
     %norm(Theta*s-Vel_coarse',2) < max_vel/10  ; %maybe more relaxed?
   cvx_end;
    else
      a0 = ones(n,1);
      lambda = 3000;
      %options1 = optimoptions(@fminunc,'MaxIterations',10000);
      V = @(a) norm(Theta*a-Vel_coarse',1) + lambda * norm(a,1) ;
      opts = psoptimset('MaxIter',20000);
      %opts.MeshTolerance = 1e-7;
      s = patternsearch(V,a0); %global optimization
    end
        
   
end



%Vel_reconstruct =  idct(s); % reconstruct full signal with high resolution
if (Flag_transform==1)
 Vel_reconstruct =  dct(s); % reconstruct full signal with high resolution
else
 Vel_reconstruct =  fft(s); % reconstruct full signal with high resolution
end

%%%%

%Velt = fft(Vel); % Fourier transformed signal
%PSD = Velt.*conj(Velt)/n; % Power spectral density

Velt = fft(Vel); % Fourier transformed signal
PSD = Velt.*conj(Velt)/n; % Power spectral density 


freq = n/(n)*(0:n);  %create the x-axis of frequencies in Hz
L = 1:floor(n/2);  % only plot the first half of freqs

xtrecon = fft(Vel_reconstruct,n);  % computes the (fast) discrete fourier transform
PSDrecon = xtrecon.*conj(xtrecon)/n;  % Power spectrum (how much power in each freq)




%%%%%%%  plot results %%%%%%%%%%



 nfft = 2^nextpow2(length(Vel));
 Fs =  1/ (t(2) - t(1)); %sampling frquency  %t = 0:1/Fs:T; 
 Pxx = abs(fft(Vel,nfft)).^2/length(Vel)/Fs; 
 Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
 figure;
 plot(Hpsd);
 title('Original PSD');




 nfft = 2^nextpow2(length(Vel_reconstruct));
 Fs =  1/ (t(2) - t(1)); %sampling frquency  %t = 0:1/Fs:T; 
 Pxx = abs(fft(Vel_reconstruct,nfft)).^2/length(Vel_reconstruct)/Fs; 
 Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
 figure;
 plot(Hpsd);
 title('Reconstructed PSD');





if (Flag_remove_DC)
    Vel = Vel + Mean_vel ;
    Vel_reconstruct = Vel_reconstruct +  Mean_vel ;
    Vel_coarse = Vel_coarse + Mean_vel;
end


 
figure;
plot(t,Vel,'LineWidth',3,'color','black');
title('Original signal','FontSize', 24);
xlabel('Time (s)','FontSize', 24);
ylabel('Flow rate cc/s','FontSize', 24);
set(gca,'fontsize',20)
set(gcf,'Position',[100 100 1000 300]);
ylim([-30 100]);
%xlim([0 0.02]);
%xlim([0.25 0.31]);

figure;
t_mat = [t(samples);Vel_coarse]';
[~,idx] = sort(t_mat(:,1));
t_mat = t_mat(idx,:);
plot(t_mat(:,1),t_mat(:,2),'LineWidth',3,'color','black');
title('Downsampled signal','FontSize', 24);
xlabel('Time (s)','FontSize', 24);
ylabel('Flow rate cc/s','FontSize', 24);
set(gca,'fontsize',20)
set(gcf,'Position',[100 100 1000 300]);
ylim([-30 100]);
%xlim([0 0.02]);

figure;
plot(t,Vel_reconstruct,'LineWidth',3,'color','black');
title('Reconsructed signal','FontSize', 24);
xlabel('Time (s)','FontSize', 24);
ylabel('Flow rate cc/s','FontSize', 24);
set(gca,'fontsize',20)
set(gcf,'Position',[100 100 1000 300]);
ylim([-30 100]);
%xlim([0 0.02]);
%xlim([0.25 0.31]);

%Gaussian filtering
wg = gausswin(5); %the higher, the more filtering
wg = wg/sum(wg);
y = filter(wg,1,Vel_reconstruct);
figure; 
plot(t,y,'LineWidth',3,'color','black');
title('Reconsructed signal with Gaussian smoothing','FontSize', 24);
xlabel('Time (s)','FontSize', 24);
ylabel('Flow rate cc/s','FontSize', 24);
set(gca,'fontsize',20)
set(gcf,'Position',[100 100 1000 300]);
ylim([-30 100]);

Error = abs( Vel_reconstruct - Vel' );
fprintf('Error in CS = %f\n', norm(Error) );

Error = abs( y - Vel' );
fprintf('Error in CS with Gaussian smoothing = %f\n', norm(Error) );

%Filter out low energy PSD data?
if(1)
     Frac_filter = 1e-2; %filter out lower than this much of max; test small values 1e-3 works!!
     xtrecon_filter = xtrecon;
     xtrecon_filter ( PSDrecon < Frac_filter*max( PSDrecon)) = 0;
     inv_xtrecon_filter = ifft(xtrecon_filter);
     figure;
     plot(t,inv_xtrecon_filter);
     title('Reconsructed signal filtered with PSD');
     %xlim([0 0.02]);    
end




