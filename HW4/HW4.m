clear all 
close all
clc

v_inf = 1;
alpha = 15;
alpha = alpha*pi/180;
s_x = -0.102;
s_y = 0;
s = s_x + 1i*s_y;
r = sqrt((1-s_x)^2+s_y^2);

% Create points for airfoil
alpha_vec = linspace(0,2*pi,1000);

for idx = 1:length(alpha_vec)
    zeta(idx) = s + r*exp(1i*alpha_vec(idx));
end

% Apply Joukowsky Transform to obtain airfoil shape
z = zeta + 1./zeta;

% Compute velocity around airfoil
gamma = 4*pi*v_inf*r*sin(alpha + asin(s_y/r));
w_tilde = v_inf.*exp(-1i.*alpha) + 1i.*gamma./(2.*pi.*(zeta - s)) - v_inf.*(r.^2)*exp(1i.*alpha)./((zeta - s).^2);
w = w_tilde./(1-(1./zeta.^2));
u = real(w);
v = -imag(w);
vel = sqrt(u.^2 + v.^2);

% Plot Airfoil
plot(real(z), -imag(z))
axis equal
l_c = max(real(z)) - min(real(z));
thickness = 100*(max(imag(z)) - min(imag(z)))/l_c;
hold on

z_t = z(1:500);
z_b = z(501:1000);
vel_t = vel(1:500);
vel_b = vel(501:1000);


% Plot pressure distribution
for j = 1:500
    cp_t(j) = 1 - (vel_t(j)./v_inf).^2;
    cp_b(j) = 1 - (vel_b(j)./v_inf).^2;
end
plot(real(z_t), -cp_t, 'm')
plot(real(z_b), -cp_b, 'c')
title('12% airfoil, 2% camber, 15 degree AoA')
xlabel('x')
ylabel('Cp')
legend('Airfoil Geometry', 'Top', 'Bottom')

% Plot coeff of lift vs. AoA
figure(2)
aoa = [0 5 10 15];
C_L_nocamb = [0 0.58 1.19 1.78];
C_L_camb = [0.47 1.064 1.65 2.23];
C_L_approx = [0 0.55 1.097 1.64];
plot(aoa, C_L_nocamb, 'r')
hold on
plot(aoa, C_L_approx, 'b')
legend('C_L', 'C_L approx')
xlabel('alpha (degrees)')
ylabel('C_L')
title('My Airfoil (2% Camber)')


