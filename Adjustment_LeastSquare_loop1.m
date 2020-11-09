%Least Square Adjustment for Traverse with Triangle Network
%Data Corrected by Party4
%This code is work for Loop1
clc
syms E38 N38 E28 N28 Ea Na
Ecu = 665580.855; %Known coordinate point
Ncu = 1519282.156;

%Observation Equation
%atan must be careful about angle
%tan( <90 degree )
%Then draw the map before made an equation.
% Azimuth CU09-05 is 313.3831267 Degree.
v1 = sqrt((Ecu-E38)^2 + (Ncu-N38)^2)-126.804;
v2 = sqrt((E38-E28)^2 + (N38-N28)^2)-56.184;
v3 = sqrt((Ecu-E28)^2 + (Ncu-N28)^2)-105.939;
v4 = atan((N38-Ncu) / (Ecu-E38))- (313.3831267-270)*pi/180.0;
v5 = atan((N28-Ncu) / (Ecu-E28)) -(107.3746776-90)*pi/180.0;
v6 = atan((E38-E28) / (N38-N28)) -(189.1653026-180)*pi/180.0;

fL = [v1;v2;v3;v4;v5;v6];
%Start value for put in adjustment loop
%These value from compass rule adjustment , It's quite accurate.
iE38 = 665488.6977 ;
iN38 = 1519369.26 ;
iE28 = 665479.749 ;
iN28 = 1519313.794 ;

mA  = jacobian(fL,[E38 N38 E28 N28]);
%Weight matrix made from SD from surveying 3 sets 
%first three are Distance SD 
%last three are Angle SD
W = diag([1.0/(0.000745355992503489)^2 1.0/(0.000687184270937547)^2 1.0/(0.000687184270939558)^2 1/(0.00028539*pi/180)^2 1/(0.000267959*pi/180)^2 1/(0.000130946*pi/180)^2]);

for i = 1:8
    A = double(subs(mA,[E38 N38 E28 N28],[iE38, iN38, iE28, iN28]));
    L = double(subs(fL,[E38 N38 E28 N28],[iE38, iN38, iE28, iN28]));
    %subs = replace parameter by value
    N = A'*W*A;
    U = A'*W*L;
    dx = -1*N\U;
    iE38 = iE38 + dx(1);
    iN38 = iN38 + dx(2);
    iE28 = iE28 + dx(3);
    iN28 = iN28 + dx(4);
    % This line for checking convergent or not
    %fprintf("E38=%.3f, N38=%.3f E28=%.3f N28=%.3f\n",dx(1),dx(2),dx(3),dx(4))
end
fprintf("Coordinate in WGS84/UTM47N (m) \n")
fprintf("E38 = %.3f, N38 = %.3f E28 = %.3f N28 = %.3f\n",iE38,iN38,iE28,iN28);
var_cov = diag(inv(N));
std = sqrt(var_cov);
fprintf("STD-E38 = %.4f, STD_N38 = %.4f STD_E28 = %.4f STD_N28 = %.4f\n",std(1),std(1),std(2),std(3));
