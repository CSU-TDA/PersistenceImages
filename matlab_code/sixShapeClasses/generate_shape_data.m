function [ToyData]=generate_shape_data(n,m,V)
% Generates 6 classes of toy data to use in classification experiments. See
% the paper "Persistence Images: A Stable Vector Representation of
% Persistent Homology"

% n is the number of points in each data set
% m is the number of examples of each type
% V is the variance of noise

ToyData=cell(m+1,6);
% Column 1: Random Point cloud (in R^3)
% Column 2: circle
% Column 3: sphere
% Column 4: 3 clusters
% Column 5: 3 clusters each with 3 smaller clusters
% Column 6: torus

% Random points drawn from a uniform distribution
ToyData{1,1}=['Random Cloud'];
for i=2:m+1
   X=rand(3,n);
   ToyData{i,1}=X;  
end

% Generate point cloud data randomly sampled from a circle in R^2 centered
% at (.5,.5) radius .5
ToyData{1,2}=['Circle'];
for i=2:m+1
   X=randn(2,n); % Guassian distribution gives even random distrubution over sphere
   s2=sum(X.^2,1);
   S=repmat(sqrt(s2),2,1);
   X=X./(2*S); % normamlize so that the vectors are all length 1/2
   X=X+repmat(.5,2,n); % shift sphere from origin to (.5,.5)
   Noise=V*randn(2,n); % generate noise, normal distribution, variance V
   X=X+Noise; 
   %scatter(X(1,:), X(2,:))
   ToyData{i,2}=X;  
end

% Generate point cloud randomly sampled from a sphere, radius 1/2, centered
% at (.5,.5,.5) with noise added
ToyData{1,3}=['Sphere'];
for i=2:m+1    
    X=randn(3,n); % Guassian distribution gives even random distrubution over sphere
    s2=sum(X.^2,1);
    S=repmat(sqrt(s2),3,1);
    X=X./(2*S); % normamlize so that the vectors are all length 1/2
    X=X+repmat(.5,3,n); % shift sphere from origin to (.5,.5,.5)
    Noise=V*randn(3,n); % generate noise, normal distribution, variance V
    X=X+Noise; 
    %scatter3(X(1,:), X(2,:),X(3,:)); 
    ToyData{i,3}=X; 
end


% clusters of points around three randomly chosen points in [0,1]x[0,1]x[0,1]
% Note: random number of points in each cluster, at most half of all points
% in a single cluster, may not be cleanly separated

% 3 tight clusters of points in 3 larger clusters around same centers as first set of clusters 

ToyData{1,4}=['Clusters'];
ToyData{1,5}=['Clusters within Clusters'];

for i=2:m+1
    Centers=rand(3);  
    a=max(randi_helper(floor(.45*n)),floor(.1*n)); % random number of points for first cluster, can be at most 45% of all points, but is at least 10%
    b=max(randi_helper(floor(.45*(n-a))), floor( (n-a)*.1)); % random number of pts in 2nd cluster, at least 10% remaining, at most 45%
    c=n-a-b; 
    X1=cat(2,repmat(Centers(:,1),1,a),repmat(Centers(:,2),1,b),repmat(Centers(:,3),1,c));  
    Y=.05*randn(3,n);
    X1=X1+Y; 
    %scatter3(X1(1,:), X1(2,:),X1(3,:)); 
    ToyData{i,4}=X1; 
 
    Centers=cat(2,repmat(Centers(:,1),1,3),repmat(Centers(:,2),1,3),repmat(Centers(:,3),1,3));  
    centers=.05*randn(3,9); %small cluster centers: perturbation from large cluster center
    Centers=Centers+centers; %center of cluster + perturbation
 
    % random number of points for each small cluster, can be at most 45% and at least 10%
    a1=max(randi_helper(floor(.45*a)),floor(.1*a)); 
    a2=max(randi_helper(floor(.45*(a-a1))),floor((a-a1)*.1)); 
    a3=a-a1-a2;    
    
    X=cat(2,repmat(Centers(:,1),1,a1),repmat(Centers(:,2),1,a2),repmat(Centers(:,3),1,a3));
       
    b1=max(randi_helper(floor(.45*b)),floor(b*.1)); % random number of points for first small cluster, can be at most 45% of all points and at least 10%
    b2=max(randi_helper(floor(.45*(b-b1))),floor((b-b1)*.1)); 
    b3=b-b1-b2;  
    
    Y=cat(2,repmat(Centers(:,4),1,b1),repmat(Centers(:,5),1,b2),repmat(Centers(:,6),1,b3));
    
    c1=max(randi_helper(floor(.45*c)),floor(c*.1)); %random number of points for first small cluster, can be at most 45% of all points and at least 10%
    c2=max(randi_helper(floor(.45*(c-c1))), floor((c-c1)*.1)); 
    c3=c-c1-c2;  
    
    Z=cat(2,repmat(Centers(:,7),1,c1),repmat(Centers(:,8),1,c2),repmat(Centers(:,9),1,c3));
    
    X=cat(2,X,Y,Z);  
    N=.02*randn(3,length(X)); %randomly disperse from centers
    X=X+N; 
    %scatter3(X(1,:), X(2,:),X(3,:)); 
    ToyData{i,5}=X;     
    
end

% Generate point cloud data randomly sampled from a torus in R^3 cantered
% at (.5,.5,.5) R=.5, r=.25
ToyData{1,6}=['Torus'];
for i=2:m+1
   [theta, phi]=torusreject(n); 
   X=[(.35+.15.*cos(theta)).*cos(phi), (.35+.15.*cos(theta)).*sin(phi),sin(theta)];
   X=X'+repmat(.5,3,n); % shift sphere from origin to (.5,.5,.5)
   Noise=V*randn(3,n); % generate noise, normal distribution, variance V
   X=X+Noise; 
   %scatter3(X(1,:), X(2,:),X(3,:)); 
   ToyData{i,6}=X;  
end


end
%save('toydata_n1','toydata_n1');
