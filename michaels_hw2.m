%{

Sources at (3,10), (9,6), and (3,2)

V_jet = 40 ft/s

A_jet = 1.77 ft^2

M_dot_jet = 70.8 ft^3/s

A_Cylinder = 7.07 ft^2

V_Cylinder = 10.0 ft/s

lamda = 295.9

V_r = lamda/(2*pi*r)

V_r = 47.1/r

V_theta = 0

 

Need to convert to cartesian coordinates

 

Source Coordinates (x1,y1)

V_r = 47.1/sqrt((x-x1)^2 + (y-y1)^2)

%}

 

 

 

% V_r Magnitude

 

        % For the Source at (3,10)

       

            % Set up the source coordinates

            x1 = 4;
            y1 = 11;

       

            % Set up the matrix

            M1 = zeros(13,16);

       

            % Use for loop to fill V_r matrix

            for x = 1:16

                    for y = 1:13
                        M1(y,x) = 47.1/sqrt((x-x1)^2 + (y-y1)^2);

                    end

            end

       

        % For the Source at (9,6)

       

            % Set up the source coordinates

            x2 = 10;

            y2 = 7;

           

            % Set up the matrix

            M2 = zeros(13,16);

           

            % Use for loop to fill V_r matrix

            for x = 1:16

                    for y = 1:13

                        M2(y,x) = 47.1/sqrt((x-x2)^2 + (y-y2)^2);

                    end

            end

       

        % For the Source at (3,2)

       

            % Set up the source coordinates

            x3 = 4;

            y3 = 3;

           

            % Set up the matrix

            M3 = zeros(13,16);

           

            % Use for loop to fill V_r matrix

            for x = 1:16

                    for y = 1:13

                        M3(y,x) = 47.1/sqrt((x-x3)^2 + (y-y3)^2);

                    end

            end

        

%          M4 = M1 + M2 + M3;

%          imagesc(M4)

 

% V_r Unit Vector X-Direction

 

        % For the Source at (3,10)

        M5 = zeros(13,16);

        for x = 1:16

                    for y = 1:13

                        M5(y,x) = (x-x1)/sqrt((x-x1)^2 + (y1-y)^2);

                    end

        end

 

        % For the Source at (9,6)

        M6 = zeros(13,16);

        for x = 1:16

                    for y = 1:13

                        M6(y,x) = (x-x2)/sqrt((x-x2)^2 + (y2-y)^2);

                    end

        end

 

        % For the Source at (3,2)

        M7 = zeros(13,16);

        for x = 1:16

                    for y = 1:13

                        M7(y,x) = (x-x3)/sqrt((x-x3)^2 + (y3-y)^2);

                    end

        end

 

% V_r Unit Vector Y-Direction

 

        % For the Source at (3,10)

        M8 = zeros(13,16);

        for x = 1:16

                    for y = 1:13

                        M8(y,x) = (y1-y)/sqrt((x-x1)^2 + (y1-y)^2);

                    end

        end

 

        % For the Source at (9,6)

        M9 = zeros(13,16);

        for x = 1:16

                    for y = 1:13

                        M9(y,x) = (y2-y)/sqrt((x-x2)^2 + (y2-y)^2);

                    end

        end

 

        % For the Source at (3,2)

        M10 = zeros(13,16);

        for x = 1:16

                    for y = 1:13

                        M10(y,x) = (y3-y)/sqrt((x-x3)^2 + (y3-y)^2);

                    end

        end

 

        % Next step is to multiply the unit vectors by the V_r magnitudes

       % for each source

 

        % After that, I must add all the x matrices together and add all

        % the y matrices together to get the u and v matrices

 

        % Last, I can quiver plot the result

 

% Magnitude adjusted Unit Vectors
       

        % For the Source at (3,10)

        U1 = M1.*M5;

        V1 = M1.*M8;

 

        % For the Source at (9,6)

        U2 = M2.*M6;

        V2 = M2.*M9;

       

        % For the Source at (3,2)

        U3 = M3.*M7;

        V3 = M3.*M10;

       

        U = U1 + U2 + U3;

        V = V1 + V2 + V3;

 


% Quiver plot

 

[Y,X]=meshgrid(1:size(U1,2),1:size(U1,1));

 

tiledlayout(1,3)

 

nexttile

 

% The Dots

plot(Y,X,'k.','MarkerSize',10)

pbaspect([1 1 1])

set(gcf, 'Position',  [50, 100, 1800, 800])

hold on

 

% The Arrows

quiver(Y,X,U,-V,'r')

title('Zero Wind Condition')

 

% The Ship

P = [1,1;1,13;16,7];
T = delaunayTriangulation(P);
triplot(T,'b')

 

 

% Now the next step is to add 10 knots of wind from the nose and then 10

% knots of wind from the starboard side.

 

% 10 knots = 16.9 ft/s

 

% The goal is going to be to construct A U matrix for the head wind and a V

% matrix for the wind on the starboard side

 

% The unit vectors are all going to be the same because the wind is uniform

% this means I can construct a vector of ones and then just multiply it by

% the magnitude and direction to get the U and V vectors.

 

U4 = -16.9*ones(13,16);

V4 = 16.9*ones(13,16);

 

% Headwind

 

U_Headwind = U1 + U2 + U3 + U4;

 

nexttile

 

% The Dots

plot(Y,X,'k.','MarkerSize',10)

pbaspect([1 1 1])

hold on

 

% The Arrows

quiver(Y,X,U_Headwind,-V,'r')

title('10 Knots from the Nose')

 

% The Ship

P = [1,1;1,13;16,7];

T = delaunayTriangulation(P);

triplot(T,'b')

 

 

% Starboard Wind

 

V_Starboard = V1 + V2 + V3 + V4;

 

nexttile

 

% The Dots

plot(Y,X,'k.','MarkerSize',10)

pbaspect([1 1 1])

hold on

 

% The Arrows

quiver(Y,X,U_Headwind,-V_Starboard,'r')

title('10 Knots from the Starboard Side')

 

% The Ship

P = [1,1;1,13;16,7];

T = delaunayTriangulation(P);

triplot(T,'b')

% Part B
% Find average over (U,-V), (U,-V_Starboard), and (U_Headwind,-V))
A_N = sqrt((mean(U,'all','omitnan'))^2 + (mean(-V,'all','omitnan'))^2);
% disp((mean(U,'all','omitnan'))^2)
disp('halu')
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