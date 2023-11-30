function [cl360, cd360] = ExtrapolationFn(alpha_data, cl, Ap, Bp, Am, Bm, thickness, camber, yc)

    cl360 = zeros(181, 1);
    cd360 = zeros(181, 1);
    
    tmpcl = [cl(dsearchn(cl, 0)-1); cl(dsearchn(cl, 0)); cl(dsearchn(cl, 0)+1)];
    tmpalpha = [alpha_data(dsearchn(cl, 0)-1); alpha_data(dsearchn(cl, 0)); alpha_data(dsearchn(cl, 0)+1)];
    alphazero = interp1(tmpcl, tmpalpha, 0);                                % alpha at CL=0
    CLzero = interp1(alpha_data, cl, 0);                                    % Cl at alpha=0
    
    if min(alpha_data) < 0 && max(alpha_data) > 0                           % Compute slope of potential graph
        slope = (cl(dsearchn(alpha_data, 0) + 1) - cl(dsearchn(alpha_data, 0) - 1))/(alpha_data(dsearchn(alpha_data, 0) + 1) - alpha_data(dsearchn(alpha_data, 0) - 1));
    elseif min(alpha_data) >= 0
        slope = (cl(2) - cl(1))/(alpha_data(2) - alpha_data(1));
    elseif max(alpha_data) <= 0
        slope = (cl(end) - cl(end-1))/(alpha_data(end) - alpha_data(end-1));
    end
    
    %% Positive Polar Extrapolation
    %{
    index_clmax = dsearchn(cl, max(cl));
    CLmax = max(cl);
    deltaalpha = alpha_data(2) - alpha_data(1);
    
    if index_clmax + Ap + 1 < size(alpha_data,1) % Ap must be integer
        a1plus = alpha_data(index_clmax + Ap);
        CL1plus = cl(index_clmax + Ap);
    else
        a1plus = (index_clmax + Ap)*deltaalpha;
        CL1plus = PlateFlow(alphazero, CLzero, a1plus, thickness, camber, yc) + 0.03;
    end
    
    if size(alpha_data,1) - (index_clmax + Bp + Ap) > 0 % Bp must be integer
        a2plus = alpha_data(index_clmax + Bp + Ap);
        CL2plus = cl(index_clmax + Bp + Ap);
    else
        a2plus = (index_clmax + Bp + Ap)*deltaalpha;
        CL2plus = PlateFlow(alphazero, CLzero, a2plus, thickness, camber, yc) + 0.03;
    end
    %}
    
    a1plus = (Ap/15+CLzero)/slope + 4;
    CL1plus = Ap/15;
    a2plus = a1plus + Bp*2;
    CL2plus = PlateFlow(alphazero, CLzero, a2plus, thickness, camber, yc) + 0.03;
    
    f1plus = ((CL1plus - PlateFlow(alphazero, CLzero, a1plus, thickness, camber, yc))/(PotFlow(CLzero, slope, a1plus) - PlateFlow(alphazero, CLzero, a1plus, thickness, camber, yc)));
    f2plus = ((CL2plus - PlateFlow(alphazero, CLzero, a2plus, thickness, camber, yc))/(PotFlow(CLzero, slope, a2plus) - PlateFlow(alphazero, CLzero, a2plus, thickness, camber, yc)));
    
    G = (abs((1/f1plus - 1)/(1/f2plus - 1)))^(1/4);
    alpha_m = (a1plus - G*a2plus)/(1 - G);
    k = (1/f2plus - 1)/(a2plus - alpha_m)^4;
    
    alpha = -1;
    for i = 1:91
        if alpha < 45
            alpha = alpha + 1;
        else
            alpha = alpha + 3;
        end
    
        if alpha < alpha_m
            delta = 0;
        else
            delta = alpha_m - alpha;
        end
    
        f = 1/(1 + k*(delta^4));
        deltaCD = 0.13*((f - 1)*PotFlow(CLzero, slope, alpha) - (1 - f)*PlateFlow(alphazero, CLzero, alpha, thickness, camber, yc));
    
        if deltaCD <= 0 
            deltaCD = 0;
        end
        
        cl360(90 + i) = f*PotFlow(CLzero, slope, alpha) + (1 - f)*PlateFlow(alphazero, CLzero, alpha, thickness, camber, yc);
        cd360(90 + i) = f*(deltaCD + 0.006 + 1.25*(thickness^2)/180*abs(alpha)) + (1 - f)*CDPlate(alpha, thickness, camber, yc) + 0.006;
    end
    
    %% Negative Polar Extrapolation
    
    a1minus = (-Am/15-CLzero)/slope - 4;
    CL1minus = -Am/15;
    a2minus = a1minus - Bm*2;
    CL2minus = PlateFlow(alphazero, CLzero, a2minus, thickness, camber, yc) - 0.03;
    
    f1minus = ((CL1minus - PlateFlow(alphazero, CLzero, a1minus, thickness, camber, yc))/(PotFlow(CLzero, slope, a1minus) - PlateFlow(alphazero, CLzero, a1minus, thickness, camber, yc)));
    f2minus = ((CL2minus - PlateFlow(alphazero, CLzero, a2minus, thickness, camber, yc))/(PotFlow(CLzero, slope, a2minus) - PlateFlow(alphazero, CLzero, a2minus, thickness, camber, yc)));
    
    G = (abs((1/f1minus - 1)/(1/f2minus - 1)))^(1/4);
    alpha_m = (a1minus - G*a2minus)/(1 - G);
    k = (1/f2minus - 1)/(a2minus - alpha_m)^4;
    
    alpha = 0;
    for i = 1:90
        if alpha > -45
            alpha = alpha - 1;
        else
            alpha = alpha - 3;
        end
    
        if alpha > alpha_m
            delta = 0;
        else
            delta = alpha_m - alpha;
        end

        f = 1/(1 + k*(delta^4));
        deltaCD = 0.13*((f - 1)*PotFlow(CLzero, slope, alpha) - (1 - f)*PlateFlow(alphazero, CLzero, alpha, thickness, camber, yc));
        
        if deltaCD <= 0
            deltaCD = 0;
        end
        
        cl360(91 - i) = f*PotFlow(CLzero, slope, alpha) + (1 - f)*PlateFlow(alphazero, CLzero, alpha, thickness, camber, yc);
        cd360(91 - i) = f*(deltaCD + 0.006 + 1.25*(thickness^2)/180*abs(alpha)) + (1 - f)*CDPlate(alpha, thickness, camber, yc) + 0.006;
    end

end

function res = PotFlow(CLzero, slope, alpha)                                % Potential flow graph(t) function 
    res = CLzero + slope*alpha;
end

function res = PlateFlow(alphazero, CLzero, alpha, thickness, camber, yc)   % Plate flow graph(s) function 
    res = (1 + CLzero/sin(pi/4)*sin(alpha*pi/180)) * CD90(alpha, thickness, camber, yc) ...
        * sin((alpha - 57.6*0.08*sin(alpha*pi/180) - alphazero*cos(alpha*pi/180))*pi/180) ...
        * cos((alpha - 57.6*0.08*sin(alpha*pi/180) - alphazero*cos(alpha*pi/180))*pi/180);
end

function res = CDPlate(alpha, thickness, camber, yc)                        % Plate CD graph function
    res = CD90(alpha, thickness, camber, yc)*(sin(alpha*pi/180))^2;
end

function res = CD90(alpha, thickness, camber, yc)                           % CD at alpha = 90 deg function
    res = 2.086 - 4.6313*yc - 1.46*thickness/2 + 1.46*camber*sin(alpha*pi/180);
end
