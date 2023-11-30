function [T, Q, P] = BEMFn(data_alpha, cl1, cl2, cd1, cd2, Block_data, HubR, PropR, BladeN, WS, RPM, Pitch_deg, ElementN, epsilon, relax, iter_limit, rho, boolPL)
    
    r_r = Block_data(1,1) + HubR;                                           % Inner section r position
    r_t = Block_data(2,1) + HubR;                                           % Outer section r position
    c_r = Block_data(1,2);                                                  % Inner section chord length
    c_t = Block_data(2,2);                                                  % Outer section chord length
    theta_r_deg = Block_data(1,3);                                          % Inner section pitch angle (deg)
    theta_t_deg = Block_data(2,3);                                          % Outer section pitch angle (deg)
    
    c = DivEle(c_r, c_t, ElementN);
    ratio = DivEle(0, 1, ElementN);
    
    dT = zeros(1, ElementN);
    dQ = zeros(1, ElementN);
    a_a = zeros(1, ElementN);
    a_t = zeros(1, ElementN);
    
    %% Compute initial aoa of each element
    theta = (DivEle(theta_r_deg, theta_t_deg, ElementN) + Pitch_deg) * (2*pi)/360;  % Pitch
    r = DivEle(r_r, r_t, ElementN);
    dr = abs(r(2) - r(1));
    
    angularspeed = (RPM/60)*2*pi;
    %lambda = PropR*angularspeed/WS;                                        % Tip Speed Ratio
    lambda_r = r*angularspeed/WS;                                           % Local Tip Speed Ratio
    
    alpha = zeros(1, ElementN);                                             % Angle of attack
    phi = zeros(1, ElementN);                                               % Wind direction
    for n = 1:ElementN
        phi(n) = atan(1/lambda_r(n));
        alpha(n) = theta(n) - phi(n);
    end
    
    %% Calculate dT(Trust), dQ(Torque) of each element
    for n = 1:ElementN
        sigma = c(n)*BladeN/(2*pi*r(n));
        
        % Interpolate cl & cd data of both section
        data_cl = cl1*(1-ratio(n)) + cl2*ratio(n);
        data_cd = cd1*(1-ratio(n)) + cd2*ratio(n);
    
        % Compute axial and tangential force coefficients
        for i = 1:iter_limit
            cl = interp1(data_alpha, data_cl, alpha(n)*180/pi);
            cd = interp1(data_alpha, data_cd, alpha(n)*180/pi);
    
            Ca = cl*cos(phi(n)) - cd*sin(phi(n));
            Ct = cl*sin(phi(n)) + cd*cos(phi(n));
    
            CN = sigma*Ca*((1 + a_a(n))/(sin(phi(n))))^2;
            
            % Prandtl Loss Correction
            F = 1;
            if boolPL
                gt = (r(n)/PropR)*tan(phi(n));
                ft = (BladeN/2)*(1 - r(n)/PropR)/gt;
                Ft = (2/pi)*acos(exp(-ft));
                
                gr = (HubR/r(n))*tan(phi(n));
                fr = (BladeN/2)*(1 + HubR/r(n))/gr;
                Fr = (2/pi)*acos(exp(-fr));
    
                F = F*Ft*Fr;
            end
    
            % Calculate induction factors
            if CN <= 0.96*F
                new_a_a = 1/(4*F*(sin(phi(n))^2)/(sigma*Ca) - 1);
            else
                new_a_a = (18*F-20-3*sqrt(abs(CN*(50-36*F)+12*F*(3*F-4))))/(36*F-50);
            end
            new_a_t = 1/(4*(sin(phi(n))*cos(phi(n)))/(sigma*Ct) + 1);
            
            if i < 11
                new_a_a = a_a(n) + 0.2*(new_a_a - a_a(n));
            else
                new_a_a = a_a(n) + relax*(new_a_a - a_a(n));
            end
    
            % Define new aoa
            phi(n) = atan((1+new_a_a)/(lambda_r(n)*(1-new_a_t)));
            alpha(n) = theta(n) - phi(n);
            
            if abs(new_a_a - a_a(n)) < epsilon && abs(new_a_t - a_t(n)) < epsilon
                break;
            end
    
            a_a(n) = new_a_a;
            a_t(n) = new_a_t;
        end
    
        Ca = cl*cos(phi(n)) - cd*sin(phi(n));
        Ct = cl*sin(phi(n)) + cd*cos(phi(n));
    
        dT(n) = Ca*rho*c(n)*(WS*(1 + a_a(n))/sin(phi(n)))^2*dr/2;
        dQ(n) = Ct*rho*c(n)*(WS*(1 + a_a(n))/sin(phi(n)))^2*r(n)*dr/2;
    
    end
    
    %% Trust, Torque and Power
    T = sum(dT)*BladeN;
    Q = sum(dQ)*BladeN*2*pi;
    P = sum(dQ)*BladeN*angularspeed;

end

function res = DivEle(element_r, element_t, number)                         % Element dividing function
    temp = linspace(element_r, element_t, number + 1);
    res = zeros(1, number);
    
    for k = 1:number
        res(k) = (temp(k) + temp(k + 1))/2;
    end
end
