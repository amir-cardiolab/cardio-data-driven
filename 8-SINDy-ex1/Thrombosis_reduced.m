function dx = Thrombosis_reduced(t,x,K_param,IC)
%Papadopoulos et al. 2014 thrombosis ODE model
%Parameters: K_param =[ K_in ; K_II^AP ; K_surf; K_AP^AP; K_AP^IIa]
%States: x = [IIa; II; AP; RP] 
%reduce X3 and X4 into one variable: X3+X4 = C = X3(0) + X4(0)

%X4 = C - X3

C = IC(3) + IC(4);

dx = [
-K_param(1)*x(1) + (K_param(3)+K_param(2)*x(3) )* x(2)  ;  %Thrombin (IIa)
-(K_param(3)+K_param(2)*x(3) )* x(2)  ; %Prothrombin (II)
K_param(4)*x(3)*(C - x(3) ) +  K_param(5)*(C - x(3) ); %Activated Platelets (AP)
%-K_param(4)*x(3)*x(4) - K_param(5)*x(4)  ; %Resting Platelets (RP)
];



end

