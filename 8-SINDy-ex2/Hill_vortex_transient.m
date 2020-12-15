function RHS = Hill_vortex_transient(t,x)
% Falahatpisheh & Kheradvar 2015 paper
%R=
T=10;  

A = 0.5*sin(2*pi*t/T);

B = 0.4*sin(2*pi*t/T);

r = sqrt( x(1)^2 + (x(2)-A)^2 + (x(3))^2 );

    
%if (r<= 1)
 u = x(1)^2 + 1 - 2*r^2;
 v = x(1)*x(2);
 %w = x(1)*x(3);
 w = x(1)*(x(3) + B);
 
%else
% u = (x(3)^2)*(r^-5) - 1/3* (r^-3) - 2/3;
% v = x(1)*x(2)*(r^-5);
% w = x(1)*x(3)*(r^-5);   
    
%end

RHS = [u; v; w];



end