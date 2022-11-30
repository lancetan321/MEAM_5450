function [gamma,Cl,camber,X,Y,U,V,panel_o_g] = computePanelData(v_inf,y0,alpha_d,N,grid_res)

% initializing variables
    alpha_r = alpha_d*(pi/180);
    v_inf_x = v_inf*cos(alpha_r);
    v_inf_y = v_inf*sin(alpha_r);
    rho = 1.225;
                

% creating the airfoil (from circle)
                
    % airfoil start point
    xi = 5; 
    yi = 0;
    
    % airfoil length and thickness
    L = 10; % do not CHANGE
    %T = for later
    
    % airfoil end point
    xf = xi+L;
    yf = yi;
    
    % circle center x coordinate
    x0 = L/2;
    
    % calculate r
    r = sqrt((xi-x0)^2+(yi-y0)^2);
    
    % initial angle
    theta_i = 0.5*(180-(2*atand(x0/y0)));
    
    % final angle 
    theta_f = 90+(90-theta_i);
    
    % scroll through angles 
    theta_vect = linspace(theta_i,theta_f,N+1)';
    
    % create vector of panel origin points
    panel_o_g = [r*cosd(theta_vect)+x0+xi,r*sind(theta_vect)+y0]; 
    
    % We have to show the camber! JKJK, unless...
    camber = (max(panel_o_g(:,2)))/L;


% creating the panels
            
    PCoords_g = zeros(N,4); 
    
    for e=1:N
        PCoords_g(e, 1) = panel_o_g(e,1);
        PCoords_g(e, 2) = panel_o_g(e,2);
        PCoords_g(e, 3) = panel_o_g(e+1,1);
        PCoords_g(e, 4) = panel_o_g(e+1,2);
    end
    %PCoords_g

 
% Panel length NOTE THIS IS A CONSTANT SO YOU CAN REMOVE THE FOR LOOP LATER
% ON
    L = zeros(N,1);
    for e=1:N 
        L(e)=sqrt( (PCoords_g(e,1)-PCoords_g(e,3))^2 + (PCoords_g(e,2)-PCoords_g(e,4))^2 );
    end
    %L

% Collocation points in the Global Coordinate Frame
    x_g = zeros(N,2);
    for e=1:N
        x_g(e,1) = (PCoords_g(e,3)-PCoords_g(e,1))/L(e); %i_hat
        x_g(e,2) = (PCoords_g(e,4)-PCoords_g(e,2))/L(e); %j_hat
    end

    y_g = zeros(N,2);
    y_g = [0 1; -1 0]*transpose(x_g);
    y_g = -transpose(y_g);
                
% Collocation points in the Global Coordinate Frame
    Col_g = zeros(N,2);
    for e=1:N
        Col_g(e,1) = PCoords_g(e,1) + 0.75*L(e)*x_g(e,1);
        Col_g(e,2) = PCoords_g(e,2) + 0.75*L(e)*x_g(e,2);
    end
    
% CREATING Icm

    for e=1:N
        for f=1:N
            
            % Collocation point location
            f_Col_x_e = Col_g(f,1) - PCoords_g(e,1);
            f_Col_y_e = Col_g(f,2) - PCoords_g(e,2);
    
            % Panel Velocity (e frame)  
            v_f_col_x_e = (1/(2*pi)) * (f_Col_y_e / ( (f_Col_y_e)^2 + (f_Col_x_e-.25*L(e))^2 )); 
            v_f_col_y_e = (-1/(2*pi)) * ((f_Col_x_e-.25*L(e)) / ((f_Col_y_e)^2 + (f_Col_x_e-.25*L(e))^2 ));
            
            % Panel Velocity (g frame)
            B = inv([x_g(e,1) y_g(e,1); x_g(e,2) y_g(e,2)]);
            A = [v_f_col_x_e v_f_col_y_e];
            C = A*B;

            % Icm
            Icm(f,e) = dot(C,y_g(f,:));
        end
    end

    %Icm

 % Creating V infinity normal

    v_inf_norm = zeros(N,1);
    for g=1:N
        v_inf_norm(g) = -dot([v_inf_x, v_inf_y], y_g(g,:));
    end

    v_inf_norm;

% calculate gamma

    gamma = inv(Icm)*v_inf_norm;

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
                  x = XX(i) - panel_o_g(j,1);
                  y = YY(h) - panel_o_g(j,2);

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