function [alpha_360, cl360, cd360] = ExtrapolationSave(alpha_data, cl, cd, cl_extra, cd_extra)

alpha_ext = [-180:3:-45 -44:44 45:3:180]';

% Combine xfoil data and extrapolation polar
alpha_360 = [alpha_ext(1:dsearchn(alpha_ext, ceil(min(alpha_data)-1))); alpha_data; alpha_ext(dsearchn(alpha_ext, floor(max(alpha_data)+1)):end)];
cl360 = [cl_extra(1:dsearchn(alpha_ext, ceil(min(alpha_data)-1))); cl; cl_extra(dsearchn(alpha_ext, floor(max(alpha_data)+1)):end)];
cd360 = [cd_extra(1:dsearchn(alpha_ext, ceil(min(alpha_data)-1))); cd; cd_extra(dsearchn(alpha_ext, floor(max(alpha_data)+1)):end)];