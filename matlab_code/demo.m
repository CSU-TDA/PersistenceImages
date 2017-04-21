%% Example 0: Generate persistence images from a "grid-like" persistence diagram with a linear weighting function

HKs=cell(1,1);
HKs{1,1}=[0,0;0,1;0,2;0,3;0,4;1,1;1,2;1,3;1,4;2,2;2,3;2,4;3,3;3,4;4,4];
res=300;
sig=0.01;
weight_func=@linear_ramp;
params=[0.05,3.5];
 %use default setting for hard/soft bounds or specify type=0 or type=1
[ PIs ] = make_PIs(HKs, res, sig, weight_func, params);

% Plot a persistence diagram and its persistence image.
HK1=HKs{1,1};
plot(HK1(:,1),HK1(:,2),'o')
hold on;
plot([0,4],[0,4]) % plot the diagonal line
figure
imagesc(PIs{1,1})

%% Example 1: Generate persistence images with linear weighting function

HKs=cell(10,5);
for i=1:10
    for j=1:5
        b=rand(50,1);
        HKs{i,j}=[b, b+(2/3)*rand(50,1)];
    end
end
res=20;
sig=0.001;
weight_func=@linear_ramp;
params=[0.05,0.5];
% use default setting for hard/soft bounds or specify type=0 or type=1
[ PIs ] = make_PIs(HKs, res, sig, weight_func, params);

% Plot a persistence diagram and its persistence image.
HK1=HKs{1,1};
plot(HK1(:,1),HK1(:,2),'o')
hold on;
plot([0,1.5],[0,1.5]) % plot the diagonal line
figure
imagesc(PIs{1,1})

%% Example 2: Generate persistence images with bump weighting function
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
[ PIs ] = make_PIs(HKs, res, sig, weight_func, params);

% Plot a persistence diagram and its persistence image.
HK1=HKs{1,1};
plot(HK1(:,1),HK1(:,2),'o')
hold on;
plot([0,1.5],[0,1.5]) % plot the diagonal line
figure
imagesc(PIs{1,1})

%% Example 3: Generate PVs with bump weighting function
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

% use default setting for hard/soft bounds or specify type=0 or type=1
[ PIs ] = make_PVs(HKs, res, sig, weight_func, params);

% Plot a persistence diagram and its persistence image.
HK1=HKs{1,1};
plot(HK1(:,1),HK1(:,2),'o')
hold on;
plot([0,1.5],[0,1.5]) % plot the diagonal line
figure
imagesc(PIs{1,1})

%% Example 4: Making persistence images with the default settings
% the default settings are a linear weighting function over the interval [0,
% max_persistence], resolution=25, sigma=1/2*(pixel height), and hard
% boundaries.

HKs=cell(10,5);
for i=1:10
    for j=1:5
        b=rand(50,1);
        HKs{i,j}=[b,b+(2/3)*rand(50,1)];
    end
end
[ PIs ] = make_PIs(HKs);

% Plot a persistence diagram and its persistence image.
HK1=HKs{1,1};
plot(HK1(:,1),HK1(:,2),'o')
hold on;
plot([0,1.5],[0,1.5]) % plot the diagonal line
figure
imagesc(PIs{1,1})
