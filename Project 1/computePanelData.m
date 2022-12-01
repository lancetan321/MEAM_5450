function [gamma,Cl,camber,X,Y,U,V,panel_origin_g] = computePanelData(v_inf,y_0,alpha_d,N,grid_res)

% initializing variables
    alpha_r = alpha_d*(pi/180);
    v_inf_x = v_inf*cos(alpha_r);
    v_inf_y = v_inf*sin(alpha_r);
    rho = 1.225;
                

% creating the airfoil (from circle)
    L = 10; 
    x_i = 5; y_i = 0;
    x_f = x_i+L; y_f = y_i;
    x_0 = x_i + L/2;
    r = sqrt((x_i-x_0)^2+(y_i-y_0)^2);
    theta_i = 0.5 * (180-(2 * atand((x_0-x_i)/y_0-y_i)));
    theta_f = 90 + (90 - theta_i);
    theta = linspace(theta_i,theta_f,N+1)';
    panel_origin_g = [r * cosd(theta) + x_0, r * sind(theta)+y_0]; 
    camber = (max(panel_origin_g(:,2)))/L;


% creating the panels
    panel_coords_g = zeros(N,4); 
    
    for i=1:N
        panel_coords_g(i, 1) = panel_origin_g(i,1);
        panel_coords_g(i, 2) = panel_origin_g(i,2);
        panel_coords_g(i, 3) = panel_origin_g(i+1,1);
        panel_coords_g(i, 4) = panel_origin_g(i+1,2);
    end

 
% Panel length 
    L = zeros(N,1);
    for i=1:N 
        L(i)=sqrt((panel_coords_g(i,1)-panel_coords_g(i,3))^2 + (panel_coords_g(i,2)-panel_coords_g(i,4))^2);
    end

% Collocation points in the Global Coordinate Frame
    x_g = zeros(N,2);
    for i=1:N
        x_g(i,1) = (panel_coords_g(i,3)-panel_coords_g(i,1))/L(i); 
        x_g(i,2) = (panel_coords_g(i,4)-panel_coords_g(i,2))/L(i); 
    end

    y_g = zeros(N,2);
    y_g = [0 1; -1 0]*transpose(x_g);
    y_g = -transpose(y_g);
                
% Collocation points in the Global Coordinate Frame
    collocation_g = zeros(N,2);
    for i=1:N
        collocation_g(i,1) = panel_coords_g(i,1) + 0.75 * L(i) * x_g(i,1);
        collocation_g(i,2) = panel_coords_g(i,2) + 0.75 * L(i) * x_g(i,2);
    end
    
% CREATING Icm
    icm = zeros(N, N);
    for i=1:N
        for j=1:N
            % x and y distances
            x_dist = collocation_g(j,1) - panel_coords_g(i,1);
            y_dist = collocation_g(j,2) - panel_coords_g(i,2);
    
            % Panel Velocity (i frame)  
            u_p = (1/(2*pi)) * (y_dist / ((y_dist)^2 + (x_dist-.25*L(i))^2)); 
            v_p = (-1/(2*pi)) * ((x_dist-.25*L(i)) / ((y_dist)^2 + (x_dist-.25*L(i))^2 ));
            
            % Panel Velocity (global frame)
            B = inv([x_g(i,1) y_g(i,1); x_g(i,2) y_g(i,2)]);
            A = [u_p v_p];
            C = A * B;

            % icm
            icm(j,i) = dot(C,y_g(j,:));
        end
    end

 % Creating V infinity normal

    v_inf_norm = zeros(N,1);
    for g=1:N
        v_inf_norm(g) = -dot([v_inf_x, v_inf_y], y_g(g,:));
    end

% calculate gamma

    gamma = inv(icm) * v_inf_norm;

% visualize the flow

    % create the meshgrid
    XX = linspace(0,20,20*grid_res);
    YY = linspace(3,-3,6*grid_res);
    [X,Y]=meshgrid(XX,YY);
    %scatter(X,Y,'.','k')

    % calculate the velocity at every grid point

    dim_x = length(XX);
    dim_y = length(YY);
    U = zeros(dim_y, dim_x);
    V = zeros(dim_y, dim_x);

    for h=1:dim_y
        for i=1:dim_x
            % scroll to first dimension in the mesh grid (1,1),
            % which is at (0,3) 
            vel_x_g = 0;
            vel_y_g = 0;
            for j=1:N 
                  %scroll through each panel to find that panels contribution to the x and y velocity at (0,3)  
                  % Find the x and y distance between the first panels origin and the first point on the mesh grid
                  % notably, I'm not saving these
                  x = XX(i) - panel_origin_g(j,1);
                  y = YY(h) - panel_origin_g(j,2);

                  %compute and rotate u_p and v_p to the global coordinate frame

                  vel_x_p = (gamma(j)/(2*pi)) * (y/((y^2)+((x-.25*L(j))^2))); % plug in the formula for u_p
                  vel_y_p = (-gamma(j)/(2*pi)) * ((x-.25*L(j))/((y^2)+((x-.25*L(j))^2))); % plug in the formula for v_p

                    B = inv([x_g(j,1) y_g(j,1); x_g(j,2) y_g(j,2)]);
                    A = [vel_x_p vel_y_p];
                    C = A*B;
                    
                    vel_x_g = vel_x_g + C(1);
                    vel_y_g = vel_y_g + C(2);
            end

            U(h,i) = v_inf_x + vel_x_g;
            V(h,i) = v_inf_y + vel_y_g;
        end
    end
    Cl = 2*sum(gamma)/(v_inf*10);
end