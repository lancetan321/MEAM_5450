clear all 
close all
clc

v_inf = input('  Enter the Speed [m/s]: ');
v = v_inf/v_inf;
alpha = input('  enter angle of attack [deg]: ');
alpha = alpha*pi/180;
s_x = input('  Circle Origin, X_0 [m]: ');
s_y = input('  Circle Origin, Y_0 [m]: ');
s = s_x + 1i*s_y;
r = input('  Radius [m]: ');


rho = 1.225;

% TRANSFORMATION PARAMETER
lambda = s_x + sqrt(r^2-s_y^2);%% this is the "C" in the currie book

% CIRCULATION
%  beta = (alpha);
%  k = 2*r*v*sin(beta);
%  Gamma = k/(2*pi) %CIRCULATION

%COMPLEX ASYMPTOTIC SPEED 
w = v * exp(1i*alpha);

%TOLLERANCE
% toll = +5e-4;

% GENERATING MESH
x = meshgrid(-10:.2:10);
y = x';

% COMPLEX PLANE
z = x + 1i*y;

% Inside-circle points are Excluded!
for a = 1:length(x)
    for b = 1:length(y)
        if abs(z(a,b)-s) <=  r
            z(a,b) = NaN;
        end
    end
end



% JOUKOWSKI TRANSFORMATION, 
J = z+lambda^2./z;

%GRAPHIC - Circle and Joukowski Airfoil
angle = 0:.1:2*pi;
z_circle = r*(cos(angle)+1i*sin(angle)) + s;
z_airfoil = z_circle+lambda^2./z_circle;

%Lift from KUTTA JOUKOWSKI THEOREM
% L = v_inf*rho*Gamma;
% L_str = num2str(L);

%COEFFICIENT OF LIFT 
Cl=2*pi*(r./lambda)*sin(alpha+s_y./r);

%Circulation necessary to satisfy Kutta Condition

Kutta_C = 4*pi*v_inf*r*sin(alpha+asin(s_y/r));

%Chord Length 
c_length = max(real(z_airfoil))-min(real(z_airfoil));

%camber(curvature) 
h = 2*s_y;

%thickness t??
t = max(imag(z_airfoil))-min(imag(z_airfoil));

% complex POTENTIAL
%f =w*(z) + (v*exp(-(1i)*alpha)*r^2)/lambda + 1i*k*log(z);


% Coefficient of lift from circulation
%lift force, Y = pro*v_inf*circulation
%C_lift = Y/(1/2)*pro*v_inf^2*l
C_lift = 2*Kutta_C/(v_inf*c_length)    ;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

%Coefficient of lift from thickness, chord length and curvature approach in
%currie

Curr_cl = 2*pi*(1+0.77*(t/c_length))*sin(alpha+2*(h/c_length));

w_tilde = v_inf.*exp(-1i.*alpha) + 1i.*Kutta_C./(2*pi.*((r-s)-z)) + v_inf.*(r.^2)*exp(1i.*alpha)./(((r-s)-z).^2);
w = w_tilde./(1 + (z-(r-s)).^2/(s.^2)); 
value = abs(w.^2);

%aphi = (s.^2 - lambda^2)/s.^2;
%PLOTTING SOLUTION
figure(1)
hold on
%contour(real(z),imag(z),real(f),[-10:.5:10]);
fill(real(z_circle),imag(z_circle),'b')
axis equal
axis on
axis([-10 10 -10 10])
title(strcat('Circle'));


figure(2)
hold on
% contour(real(J),imag(J),real(f),[-10:.5:10])

fill(real(z_airfoil),imag(z_airfoil),'b')
axis equal
axis on
axis([-10 10 -10 10])
title('Transformed Airfoil');

P_airfoil = 99500 + 0.5*1.225*((v_inf.^2) - (abs(w).^2)); 
%Pressure distribution
 figure
 a = linspace(0,c_length,101);
 Cp = (P_airfoil - 99500)./(0.5*1.225*((v_inf.^2) - (abs(w).^2)));
 plot(a,Cp);
%Coefficient of pressure

% cp = 1 - 4*sin(t).^2 + 2* G / (pi*a*V_i) *sin(t) - (2* G/ (pi*a*V_i) )^2 ;

