function [y] = rect_region(interval,params)

y = 1;
if interval(1) < params(1) || interval(1) > params(2) || interval(2) < params(3) || interval(2) > params(4)
    y = 0;
    return
end


end

