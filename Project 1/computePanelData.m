function [gamma,C_l,camber,X,Y,U,V,panel_origin_g] = computePanelData(v_inf,y_0,alpha_d,N,grid_res)

    % initializing variables
    alpha_r = alpha_d*(pi/180);
    v_inf_x = v_inf*cos(alpha_r);
    v_inf_y = v_inf*sin(alpha_r);
    rho = 1.225;
%     rho = .00237689240;

    % creating the airfoil 
    L = 10; 
%     L = 1.5;
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

    % computing panel lengths
    L = zeros(N,1);
    for i=1:N 
        L(i)=sqrt((panel_coords_g(i,1)-panel_coords_g(i,3))^2 + (panel_coords_g(i,2)-panel_coords_g(i,4))^2);
    end

    % computing collocation points in global frame
    x_g = zeros(N,2);
    for i=1:N
        x_g(i,1) = (panel_coords_g(i,3)-panel_coords_g(i,1))/L(i); 
        x_g(i,2) = (panel_coords_g(i,4)-panel_coords_g(i,2))/L(i); 
    end

    y_g = zeros(N,2);
    y_g = [0 1; -1 0]*transpose(x_g);
    y_g = -transpose(y_g);
                
    collocation_g = zeros(N,2);
    for i=1:N
        collocation_g(i,1) = panel_coords_g(i,1) + 0.75 * L(i) * x_g(i,1);
        collocation_g(i,2) = panel_coords_g(i,2) + 0.75 * L(i) * x_g(i,2);
    end
    
    % creating icm matrix
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

            % updating icm matrix entry
            icm(j,i) = dot(C,y_g(j,:));
        end
    end

    % computing v_infinity normal and calculating gamma
    v_inf_norm = zeros(N,1);
    for i=1:N
        v_inf_norm(i) = -dot([v_inf_x, v_inf_y], y_g(i,:));
    end
    gamma = inv(icm) * v_inf_norm;

    % visualizing  flow
    mesh_X = linspace(0,20,20*grid_res);
    mesh_Y = linspace(3,-3,6*grid_res);
    [X,Y]=meshgrid(mesh_X,mesh_Y);
    dim_x = length(mesh_X);
    dim_y = length(mesh_Y);
    U = zeros(dim_y, dim_x);
    V = zeros(dim_y, dim_x);
    C_l = 2*sum(gamma)/(v_inf*10);

    for i=1:dim_y
        for j=1:dim_x
            u_p_g = 0;
            v_p_g = 0;
            for k=1:N 
                x = mesh_X(j) - panel_origin_g(k,1);
                y = mesh_Y(i) - panel_origin_g(k,2);
        
                u_p = (gamma(k)/(2*pi)) * (y/((y^2)+((x-.25*L(k))^2)));
                v_p = (-gamma(k)/(2*pi)) * ((x-.25*L(k))/((y^2)+((x-.25*L(k))^2)));
        
                B = inv([x_g(k,1) y_g(k,1); x_g(k,2) y_g(k,2)]);
                A = [u_p v_p];
                C = A * B;
                    
                u_p_g = u_p_g + C(1);
                v_p_g = v_p_g + C(2);
            end

            U(i,j) = v_inf_x + u_p_g;
            V(i,j) = v_inf_y + v_p_g;
        end
    end
%     lift = sum(gamma)*v_inf*rho;
%     disp(lift)
end