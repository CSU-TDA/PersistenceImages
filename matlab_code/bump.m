function [b] = bump(b_p_data, params)
% bump determines the approximate value at radius of a smooth, 
% rotationally-invariant bump function h(x) with the properties:
% 
%   1.) h(z) = 1 if |z| <= minradius
%   2.) 0 < h(z) < 1 if minradius < |z| < maxradius 
%   3.) h(z) = 0 if maxradius < |z|

% persistence length
radius = b_p_data(2);

% Parameters controlling range of non-zero, non-trivial weightings
maxradius = params(1);
minradius = params(2);

% Parameters controlling bump function near maxradius and minradius
M = [params(3),params(4)];

% Parameter to control error approximation.
resolution = params(5);

if radius > maxradius
   b = 1; 
   return;
elseif radius < minradius
   b = 0;
   return;
end


x = linspace(minradius, maxradius, resolution);
g = zeros(1,resolution);

for i=1:resolution
    g(i) = f(x(i)-minradius, M(1))*f(maxradius-x(i), M(2));
end

b = h(x,g,radius);



% ================
% Local functions
% ================
% Smoothly varying, monotonic function whose range is [0,1).
function y = f(x,M)
    if x > 0
        y = exp(-M/x);
    else
        y = 0;
    end


    
% Integral function h(x) = int_{-inf}^{radius} g(x) / int_{-inf}^{inf} g(x)
function z = h(x,g,radius)
    [~,idx] = min(abs(x-radius));
    if idx < 2
        z = 0;
    else
        z = trapz(x(1:idx),g(1:idx))/trapz(x,g);
    end


