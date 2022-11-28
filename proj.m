% Part 1: Define airfoil coordinates
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

% Plot Airfoil
plot(real(z), -imag(z))
axis equal
l_c = max(real(z)) - min(real(z));
thickness = 100*(max(imag(z)) - min(imag(z)))/l_c;
hold on




