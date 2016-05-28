%% A minimal example using the make_PIs and make_PVs functions
close all
clear all
load('circleAndTorusIntervals.mat')
% circleAndTorusIntervals_H1 is a 1x2 cell containing the H1 intervals for
% points sampled from a noisy circle and torus.
% circleAndTorusIntervals_H0 is a 1x2 cell containing the H0 intervals from
% points sampled from a noisy circle and torus.

sig=0.0001;
res=50;
[H1_PIs] = make_PIs(circleAndTorusIntervals_H1, res, sig);
figure, imagesc(H1_PIs{1})
title('Persistant Image for H1 diagram of points sampled from a noisy circle')
figure, imagesc(H1_PIs{2})
title('Persistant Image for H1 diagram of points sampled from a noisy torus')

% Since all of the H0 intervals begin at time zero, we produce a
% 1-dimensional "vector grid" on top of this persistence diagram (instead
% of a 2-dimesional "image grid").
[H0_PVs] = make_PVs(circleAndTorusIntervals_H0, res, sig);
figure, imagesc(H0_PVs{1}), axis equal
title('Persistant Vector for H0 diagram of points sampled from a noisy circle')
figure, imagesc(H0_PVs{2}), axis equal
title('Persistant Vector for H0 diagram of points sampled from a noisy torus')

%% Convert images to vectors

% Before applying machine learning techniques, it is often helpful to
% conver a persistent image in matrix format into a vector by concatenating
% columns. This is accomplished by the following command.
H1_vecs = vecs_from_PIs(H1_PIs);