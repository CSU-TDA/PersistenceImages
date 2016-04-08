%% Example 1: generate PIs with linear weighting function

HKs=cell(10,5);
for i=1:10
    for j=1:5
        b=rand(50,1);
        HKs{i,j}=[b, b+(2/3)*rand(50,1)];
    end
end
res=20;
sig=.01;
weight_func=@linear_ramp;
params=[.5,1.5];
%use default setting for hard/soft bounds or specify type=0 or type=1
[ PIs ] = make_PIs(HKs, res, sig, weight_func,params);

%% Example 2: generate PIs with bump weighting function
HKs=cell(10,5);
for i=1:10
    for j=1:5
        b=rand(50,1);
        HKs{i,j}=[b, b+(2/3)*rand(50,1)];
    end
end
res=20;
sig=.01;

maxradius = 1.5;
minradius = .5;
maxgrow = 1;
mingrow = 1;
resolution = 10000;

weight_func = @bump;
params = [maxradius, minradius, maxgrow, mingrow, resolution];

%use default setting for hard/soft bounds or specify type=0 or type=1
[ PIs ] = make_PIs(HKs, res, sig, weight_func,params);

%% Example 3: generate PVs with bump weighting function
HKs=cell(10,5);
for i=1:10
    for j=1:5
        HKs{i,j}=[zeros(50,1),(2/3)*rand(50,1)];
    end
end
res=20;
sig=.01;

maxradius = 1.5;
minradius = .5;
maxgrow = 1;
mingrow = 1;
resolution = 10000;

weight_func = @bump;
params = [maxradius, minradius, maxgrow, mingrow, resolution];

%use default setting for hard/soft bounds or specify type=0 or type=1
[ PIs ] = make_PVs(HKs, res, sig, weight_func,params);

%% Example 4: making PIs with the default settings
%the default settings are a linear weighting function over the interval [0,
%max_persistence], resolution=25, sigma=1/2*(pixel height), and hard
%boundaries.

HKs=cell(10,5);
for i=1:10
    for j=1:5
        b=rand(50,1);
        HKs{i,j}=[b,b+(2/3)*rand(50,1)];
    end
end
[ PIs ] = make_PIs(HKs);




