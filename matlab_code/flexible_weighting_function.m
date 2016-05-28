% Random persistence diagram
intervals = zeros(500,2);
intervals(:,1) = rand(500,1);
intervals(:,2) = intervals(:,1) + 2/3*abs(randn(500,1));

% =================
%  Example 1
% =================
% Specifying the function.  In this case it is a piecewise linear function. 
% which is identically 0 until params(1) and identically 1 after params(2).
weight_func = @linear_ramp; % <-- the user would input the function handle specifying their weighting function
params = [1/5,1.5]; % <-- the user would supply any parameters needed to specify their weighting function. 

% ------------------------------------------------------------------------------------
% fixed code that computes weight at each persistence interval
weights = arrayfun(@(row) weight_func(intervals(row,:), params), 1:size(intervals,1))';
% ------------------------------------------------------------------------------------

% plot the persistence diagram with points colored by weighting function
figure
scatter(intervals(:,1),intervals(:,2),[],weights,'filled')
axis equal
axis([0,1,0,max(intervals(:,2))])
colorbar


% =================
%  Example 2
% =================
% Specifying the function.  In this case it is our bump function with  
% a host of parameters specifying its form.

maxradius = 1.5;
minradius = .5;
maxgrow = 1;
mingrow = 1;
resolution = 10000;

weight_func = @bump; % <-- the user would input the function handle specifying their weighting function
params = [maxradius, minradius, maxgrow, mingrow, resolution];
% <-- the user would supply any parameters needed to specify their weighting function. 

% ------------------------------------------------------------------------------------
% fixed code that computes weight at each persistence interval
weights = arrayfun(@(row) weight_func(intervals(row,:), params), 1:size(intervals,1))';
% ------------------------------------------------------------------------------------

% plot the persistence diagram with points colored by weighting function
figure
scatter(intervals(:,1),intervals(:,2),[],weights,'filled')
axis equal
axis([0,1,0,max(intervals(:,2))])
colorbar


% =================
%  Example 3
% =================
% Compare the two methods of computing the linear weighting function
weight_func = @linear_ramp;
params = [1/5,1.5];
tic
linear_ramp_weights = arrayfun(@(row) weight_func(intervals(row,:), params), 1:size(intervals,1))';
linear_ramp_time = toc


weight_func = @bump;
params = [1.5, 1/5, 0, 0, 10000];
tic
linear_bump_function_weights = arrayfun(@(row) weight_func(intervals(row,:), params), 1:size(intervals,1))';
linear_bump_function_time = toc

max_weighting_function_difference = max(abs(linear_bump_function_weights - linear_ramp_weights))


% =================
%  Example 4
% =================
% Specifying the function.  In this case it is a function which is 0
% outside a rectangular region of the persistence diagram.


weight_func = @rect_region; % <-- the user would input the function handle specifying their weighting function
params = [0,1/2,1,1.5];% <-- the user would supply any parameters needed to specify their weighting function. 

% ------------------------------------------------------------------------------------
% fixed code that computes weight at each persistence interval
weights = arrayfun(@(row) weight_func(intervals(row,:), params), 1:size(intervals,1))';
% ------------------------------------------------------------------------------------

% plot the persistence diagram with points colored by weighting function
figure
scatter(intervals(:,1),intervals(:,2),[],weights,'filled')
axis equal
axis([0,1,0,max(intervals(:,2))])
colorbar