%Least Square Adjustment for Traverse with Triangle Network
%Data Corrected by Party3
%This code is work for loop2
clc
syms E38 N38 E28 N28 Ea Na
Ecu = 665580.855;
Ncu = 1519282.156;

%Observation Equation
%atan must be careful about angle
%tan( <90 degree )
%Then draw the map before made an equation.
% Azimuth CU09-05 is 313.3831267 Degree.
v1 = sqrt((Ecu-E38)^2 + (Ncu-N38)^2)-126.804;
v2 = sqrt((Ea-Ecu)^2 + (Na-Ncu)^2)-58.438;
v3 = sqrt((Ea-E38)^2 + (Na-N38)^2)-115.161;
v4 = atan((Ea-Ecu) /  (Na-Ncu) ) - 18.5112517*pi/180.0 ;
v5 =atan((N38-Na) /  (Ea-E38) ) -(285.9700248-270)*pi/180.0 ;
v6 =atan((N38-Ncu)  /(Ecu-E38) ) - (133.3837749 - 90)*pi/180.0 ;

fL = [v1;v2;v3;v4;v5;v6];
%Start value for put in adjustment loop
%These value from compass rule adjustment , It's quite accurate.
iE38 = 665488.6977 ;
iN38 = 1519369.26 ;
iE28 = 665479.749 ;
iN28 = 1519313.794 ;
iEa = 665599.409;
iNa = 1519337.57 ;
u = [E38 N38 Ea Na];
mA  = jacobian(fL,u);
%Weight matrix made from SD from surveying 3 sets 
%first three are Distance SD 
%last three are Angle SD
W = diag([1.0/(0.000942809)^2 1.0/(0.00057735)^2 1.0/(0.000816497)^2 1/(0.000290968*pi/180)^2 1/(0.000229155*pi/180)^2 1/(0.000360101*pi/180)^2]);

for i = 1:8
    A = double(subs(mA,u,[iE38, iN38, iEa, iNa]));
    L = double(subs(fL,u,[iE38, iN38, iEa, iNa]));
    %subs = replace parameter by value
    N = A'*W*A;
    U = A'*W*L;
    dx = -1*N\U;
    iE38 = iE38 + dx(1);
    iN38 = iN38 + dx(2);
    iEa = iEa + dx(3);
    iNa = iNa + dx(4);
    % This line for checking convergent or not
    %fprintf("E38=%.3f, N38=%.3f Ea=%.3f Na=%.3f\n",dx(1),dx(2),dx(3),dx(4))
end
fprintf("Coordinate in WGS84/UTM47N (m) \n")
fprintf("E38 = %.3f, N38 = %.3f Ea = %.3f Na = %.3f\n",iE38,iN38,iEa,iNa);
var_cov = diag(inv(N));
std = sqrt(var_cov);
fprintf("STD-E38 = %.4f, STD_N38 = %.4f STD_Ea = %.4f STD_Na = %.4f\n",std(1),std(1),std(2),std(3));
