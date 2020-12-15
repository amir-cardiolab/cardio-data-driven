function yout = poolData_control(yin,nVars,polyorder,t)
n = size(yin,1);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

if(polyorder>=2)    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

%sines and cosines
K_s = 5; %10
Omega = 2*pi* 0.05; 
 for k=1:K_s;
    yout = [yout sin(k*Omega*t) cos(k*Omega*t)];
 end
 
 if(1)
 %sines and cosines mixed with polynomials
 for k=1:K_s;
    yout = [yout yin(:,1).*sin(k*Omega*t) yin(:,1).*cos(k*Omega*t)];
    yout = [yout yin(:,2).*sin(k*Omega*t) yin(:,2).*cos(k*Omega*t)];
    yout = [yout yin(:,3).*sin(k*Omega*t) yin(:,3).*cos(k*Omega*t)];
 end
 
 end
 
 
 
    