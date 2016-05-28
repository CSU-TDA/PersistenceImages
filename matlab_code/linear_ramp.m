function [y] = linear_ramp(b_p_data, params)


min_radius = params(1);
max_radius = params(2);

x = b_p_data(:,2);

if x <= min_radius
    y = 0;
elseif x >= max_radius
    y = 1;
else
    y = (x - min_radius)/(max_radius - min_radius);
end



end

