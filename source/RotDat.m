function [res_x, res_z] = RotDat(x, z, angle_deg)

angle = -angle_deg*pi/180;                                                  % Degree to radian
res_x = zeros(length(x),1);
res_z = zeros(length(x),1);

for n = 1:length(x)                                                         % Rotate data
    res_x(n) = x(n)*cos(angle) - z(n)*sin(angle);
    res_z(n) = x(n)*sin(angle) + z(n)*cos(angle);
end