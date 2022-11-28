% Part B
% Find average over (U,-V), (U,-V_Starboard), and (U_Headwind,-V))
A_N = sqrt((mean(U,'all','omitnan'))^2 + (mean(-V,'all','omitnan'))^2);
disp((mean(U,'all','omitnan'))^2)
A_S = sqrt((mean(U,'all','omitnan'))^2 + (mean(-V_Starboard,'all','omitnan'))^2);
A_H = sqrt((mean(U_Headwind,'all','omitnan'))^2 + (mean(-V,'all','omitnan'))^2);
%Suction force is the difference between atmospheric pressure and the
%pressure under the aircraft.
%No Wind
P_ambient = 2116; %psf
rho = 0.0765;
V_avg = 16.87; %CHANGE CASE
V_avg_under = A_H; %CHANGE CASE
P_dynamic = 0.5*rho*V_avg^2;
P_dynamic_under = 0.5*rho*(V_avg_under)^2;
P_static_under = P_ambient + P_dynamic - P_dynamic_under;
Area = 90; %ft^2
F = (P_ambient-P_static_under)*Area
% No Wind (Force was 64.19 lbs)
% Starboard (Force was 67.68 lbs)
% Headwind (Force was -434.8 lbs)