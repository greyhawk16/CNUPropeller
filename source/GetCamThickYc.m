function [camber, thickness, yc] = GetCamThickYc(Af_dat)
    a1 = Af_dat;

    % LE, TE's Value
    TE = max(a1(:, 1));                                                     % x-value of TE
    LE = min(a1(:, 1));                                                     % x-value of LE
    
    TE_loc = find(a1(:, 1) == TE)';                                         % find all locations of TE
    LE_loc = find(a1(:, 1) == LE)';                                         % find all locations of LE
    
    locs = sort([LE_loc, TE_loc]);                                          % stores indexes of LE and TE. sort indexes in ascending order
    
    a1_size = size(a1);
    a1_rows = a1_size(1);
    
    if locs(end) < a1_rows                                                  % if set of indexes(locs) does NOT end with (no. of rows in table) 
        locs(end+1) = a1_rows;                                              % append (no. of rows in table) to locs
    end
    if locs(1) > 1                                                          % if locs does NOT start with 1
        locs = [1, locs];                                                   % add 1 to the start of locs
    end
    
    tmp = size(locs); 
    rows = tmp(2); 
    
    if rem(rows, 2) ~= 0                                                    % 3 points
        side_1 = a1(locs(1):locs(2), :);                                    % 1st-half
        side_2 = a1(locs(2):locs(3), :);                                    % 2nd-half

    else                                                                    % 4 points
        side_1 = a1(locs(1):locs(2), :);                                    % 1st-half
        side_2 = a1(locs(3):locs(4), :);                                    % 2nd-half

    end
    
    thickness = 0;
    camber = 0;
    
    % 1st-half -> 2nd-half
    tmp = size(side_1);
    ed_1 = tmp(1);
    
    for i=1:ed_1
        p1 = side_1(i, :);                                                  % point at 1st-half
        p1y = p1(1, 2);
        p2x = side_1(i, 1);                                                 % point to look at 2nd-half
        p2y = interp1(side_2(:, 1), side_2(:, 2), p2x);                     % linear interpolation
    
        thickness = max(thickness, abs(p1y - p2y));
        camber = max(camber, (p1y + p2y) / 2);
    end
    
    % 2nd-half -> 1st-half
    tmp = size(side_2);
    ed_2 = tmp(1);
    
    for i=1:ed_2
        p2 = side_2(i, :);                                                  % point at 2nd-half
        p2y = p2(1, 2); 
        p1x = side_2(i, 1);                                                 % point to look at 1st-half
        p1y = interp1(side_1(:, 1), side_1(:, 2), p1x);                     % linear interpolation
    
        thickness = max(thickness, abs(p1y - p2y));
        camber = max(camber, (p1y + p2y) / 2);
    end

    yc = abs(interp1(side_1(:, 1), side_1(:, 2), 0.0125)) + abs(interp1(side_2(:, 1), side_2(:, 2), 0.0125));   % y/c at x/c = 0.0125
    
end