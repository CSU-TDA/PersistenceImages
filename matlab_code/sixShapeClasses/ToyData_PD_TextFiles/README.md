# 6 shape classes
This subfolder contains persistence diagram data for 6 shape classes:
(1) A unit cube
(2) A circle of diameter one
(3) A sphere of diameter one
(4) Three clusters with centers randomly chosen in the unit cube
(5) Three clusters within three clusters
(6) A torus with a major diameter of one and a minor diameter of one half.
These shape classes are described in Section 6.1 of the paper "Persistence Images: A Stable Vector Representation of Persistent Homology"

We produce 25 point clouds of 500 randomly sampled points from each shape class. We then add a level of Gaussian noise to each point, at a noise level neta=0.1 or neta=0.05. We then have already computed the persistent homology intervals in homological dimension i=0 and i=1.

For example, the file 
ToyData_PD_n05_23_6_0.txt
corresponds to noise level neta=0.05, the 23rd point cloud randomly sampled from shape class (6), with persistent homology computed in dimension 0.
Each row of this file has two entries: the birth and death time of a 0-dimensional persistent homology interval.

By contrast, the file 
ToyData_PD_n1_21_3_1.txt
corresponds to noise level neta=1, the 21st point cloud randomly sampled from shape class (3), with persistent homology computed in dimension 1.
Each row of this file has two entries: the birth and death time of a 1-dimensional persistent homology interval.

Your task is to use machine learning to distinguish these six classes from each other. In a K-medoids clustering test, some accuracies and computation times are displayed for bottleneck distances, Wasserstein distances, persistence landscapes, and persistence images in Table 1 of "Persistence Images: A Stable Vector Representation of Persistent Homology". Do you have ideas for beating these accuracies or computation times?
