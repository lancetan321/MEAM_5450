% data pre-processing
airfoil_data = readcell('airfoil_data.txt');
airfoil_data = airfoil_data(2:end-1,:);
x_Upper = cell2mat(airfoil_data(:,5));
x_Lower = cell2mat(airfoil_data(:,7));
v_Upper = cell2mat(airfoil_data(:,6));
v_Lower = cell2mat(airfoil_data(:,8));
v_FreeStream = 10; % m/s

% initialize variables
P_amb = 101.3*10^3; % Pa
rho_amb = 1.225; % kg/m^3

% Total air pressure
P_tot = P_amb + 1/2*rho_amb*v_FreeStream^2;

% Calculate static air pressure
P_st_u = (P_tot - 1/2*rho_amb*v_Upper.^2)/1000; % kPa
P_st_l = (P_tot - 1/2*rho_amb*v_Lower.^2)/1000; % kPa

% plot pressure distribution
hold on
plot(x_Lower,P_st_l,'b')
plot(x_Upper,P_st_u,'m')
hold off
xlabel('x (non-dim)')
ylabel('Static Pressure (kPa)')
legend('Upper','Lower')
grid on