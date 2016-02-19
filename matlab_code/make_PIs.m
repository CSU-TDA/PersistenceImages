function [ PIs ] = make_PIs(interval_data, res, sig, bump_func_params, type)
%make_PIs generates the set of persistence images for the PH interval data
%stored in the cell array titled interval_data.
% INPUTS:   -interval_data: is a cell array containing the interval data
%            for all of the point clouds in consideration. Each sheet in
%            the cell array corresponds to a different Betti dimension.
%            Each interval set is assumed to be an nX2 matrix with the
%            first column corresponding to the birth times of features and
%            the second column is the death time of each feature. All death
%            times must be greater than the birth time and there must be a
%            finite death time for each feature. 
%            -res: the desired resolution of the image
%            -sig: is the desired variance of the Gaussians used to compute
%            the images
%            bump_func_params: a 1x2 matrix containing the values of the
%            bump function parameters [M1,M2]. M1 controls how quickly the
%            wieghting function leaves 0 and the M2 controls how quickly
%            the function levels off to 1.
%            -type: refers to the declaration of a soft' or 'hard' with
%            respect to the boundaries. type=1 produces hard bounds, type=2
%            produces soft boundaries on the images.
%OUTPUTS:   -PIs: The set of persistence images generated based on the
%            options specified for the provided interval data.
[ b_p_data, max_b_p_Hk, problems ] = birth_persistence_coordinates(interval_data );
%first do a check to make sure all the points (birth,persistence) points
%are viable.
if size(problems,1)>0
    error('Error: Negative Persistences Present')
elseif size(problems,1)==0;
    'All positive Persistence, continuing'
end

if nargin>5
    error('Error: too many input arguments')
elseif nargin==5
    res=res;
    sig=sig;
    bump_func_params=bump_func_params;
    type=type;
elseif nargin==4
    res=res;
    sig=sig;
    bump_func_params=bump_func_params;
    type=1; %default boundary setting is hard.
elseif nargin==3 
    res=res;
    sig=sig;
    bump_func_params=[0,0]; %default setting is a linear bump function.
    type=1; %the default boundary setting is hard
elseif nargin==2
    res=res;
    sig=.5*(max(max(max_b_p_HK(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
    bump_func_params=[0,0]; %default setting is a linear bump function
    type=1; %the default boundary setting is hard.
elseif nargin==1
    res=25; %default resolution is equal to 25 pixels.
    sig=.5*(max(max(max_b_p_HK(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
    bump_func_params=[0,0]; %default setting is a linear bump function
    type=1; %the default boundary setting is hard.
end
    
    
if type==1       
    [ data_images ] = hard_bound_PIs( b_p_data, max_b_p_Hk, bump_func_params, res,sig);
elseif type==2
    [ data_images ] = soft_bound_PIs( b_p_data, max_b_p_Hk, bump_func_params, res,sig);
end
    PIs=data_images;
    
function [ b_p_data, max_b_p_Hk, problems] = birth_persistence_coordinates(interval_data )
%birth_persitence_coordinates takes in the interval data as output by the
%duke TDA code (birth-death coordinates) and changes them into
%birth-persistence coordinates.
%INPUTS:      -interval_data: is a cell array containing the interval data
%            for all of the point clouds in consideration. Each sheet in
%            the cell array corresponds to a different Betti dimension.
%            Each interval set is assumed to be an nX2 matrix with the
%            first column corresponding to the birth times of features and
%            the second column is the death time of each feature. All death
%            times must be greater than the birth time and there must be a
%            finite death time for each feature. 
%OUTPUT:     -b_p_interval_data: This is the modified coordinate data
%            in a cell array. The sheets contain the modified Hk data.
%            -max_b_p_Hk: gives the maximal persistence and maximal
%            birth time across all point clouds for each Hk. 
%            This information is used to create the boundaries for the
%            persistence images.
%         
[m,n,o]=size(interval_data);
%allocate space
max_persistences=zeros(m,n,o);
max_birth_times=zeros(m,n,o);
birth_persistence=cell(m,n,o);
problems=[];
for k=1:o
for i=1:n
    for j=1:m
        B=interval_data{j,i,k};
        %pulls the Hk interval data for the (j,i)th point cloud.
        max_persistences(j,i,k)=max(B(:,2)-B(:,1));
        %computes the Hk persistence (death-birth) for the (j,i)th point cloud
        max_birth_times(j,i,k)=max(B(:,1));
        %determines that maximal birth time for an Hk feature for the
        %(j,i)th point cloud. We will take the maximum over all of point
        %clouds to generate non-normalized Hk PIs. 
        C=B(:,2)-B(:,1);
        birth_persistence{j,i,k}=[B(:,1) C];
        %birth-persistence coordinates for Hk
%        end
        D=find(C<=0);
        if length(D)>0
            problems=[problems; j,i,k];
        elseif length(D)==0
            problems=problems;
        end

    end
end
Hk_max_birth(k,1)=max(max(max_birth_times(:,:,k)));
%determine the maximum birth time of all Hk features across the point
%clouds
Hk_max_persistence(k,1)=max(max(max_persistences(:,:,k)));
%determine the maximum persistence of all Hk features across the point
%clouds
end
max_b_p_Hk=[Hk_max_birth, Hk_max_persistence];
b_p_data=birth_persistence;

end
function [lin_interp_vals] = approx_bump(X,min,max,varargin)
%approx_bump takes in a vector of values and approximates the bump function
%at these values. This can then be used to generate linear approximation at
%a novel point
lin_interp_vals=zeros(length(X),2);
% Check for threshold parameter
if size(varargin) == 0
    M = [1,1];
else
    M = varargin{1};
end

for i=1:length(X)
    lin_interp_vals(i,1)=X(i);
    lin_interp_vals(i,2)=bump(X(i), min, max, M);
end
    
end
function [b] = bump(radius, minradius, maxradius, varargin)
% bump2 determines the approximate value at radius of a smooth, 
% rotationally-invariant bump function h(x) with the properties:
% 
%   1.) h(z) = 1 if |z| <= minradius
%   2.) 0 < h(z) < 1 if minradius < |z| < maxradius 
%   3.) h(z) = 0 if maxradius < |z|

if radius > maxradius
   b = 0; 
   return;
end

% Check for threshold parameter
if size(varargin) == 0
    M = [1,1];
else
    M = varargin{1};
end

% Parameter to control error approximation.
resolution = 10000;

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
end

% Integral function h(x) = int_{-inf}^{radius} g(x) / int_{-inf}^{inf} g(x)
function z = h(x,g,radius)
    [~,idx] = min(abs(x-radius));
    if idx < 2
        z = 0;
    else
        z = trapz(x(1:idx),g(1:idx))/trapz(x,g);
    end
end

end
function [ data_images ] = hard_bound_PIs( b_p_data, max_b_p_Hk, bump_func_params, res,sig)
%hard_bound_PIs generates the PIs for a set of point clouds with the
%b_p_data. Hard refers to the fact that we cut the boundaries off hard at
%the maximum values. 
%   INPUTS:        -b_p_data: birth-persistence points for each of the
%                  point clouds. Each sheet corresponds to a different
%                  Betti dimension.
%                  -max_b_p_Hk: gives the maximal persistence and maximal
%                  birth time across all point clouds for each Hk. 
%                  This information is used to create the boundaries for 
%                  the persistence images.
%                  -bump_func_params: give the values for the
%                  bump function parameters which determine how fast the
%                  function leaves zero and how quickly it levels off at
%                  1. [M1,M2] where M1 controls the speed of leaving zero
%                  and M2 controls the speed of leveling off.
%                  -res: the resolution (number of pixels). we create
%                  square images with rectangular pixels.
%                  -sig: the variance of the gaussians. 
%   OUTPUT:        -data_images: the set of persistence images generated
%                  using the selected parameter values. Each sheet
%                  corresponds to the PIs generated for different Hk
%                  interval data.

[m,n,o]=size(b_p_data);
%allocate space
data_images=cell(m,n,o);
for k=1:o  
Hk_max_b=max_b_p_Hk(k,1);
Hk_max_p=max_b_p_Hk(k,2);    
%select 100 sample values of the bump function to linearly interpolate
%between to get scaling factors. For now, no matter what the range we only
%compute the value at 100 points.
X_Hk=linspace(0,Hk_max_p,100);    
%compute the bump function values at the test points. 
[lin_interp_vals_Hk] = approx_bump(X_Hk,0,Hk_max_p,[0,0]);    
sigma=[sig,sig]; %duplicate the variance for the PIs      
%set up gridding for Hk
birth_stepsize_Hk=Hk_max_b/res; %the x-width of a pixel
persistence_stepsize_Hk=Hk_max_p/res; %the y-height of a pixel
grid_values1_Hk=0:birth_stepsize_Hk:Hk_max_b; %must be increasing from zero to max_dist
grid_values2_Hk=Hk_max_p:-persistence_stepsize_Hk:0; %must be decreasing from max_dist to zero
            for p=1:m
                for t=1:n
                Hk=b_p_data{p,t,k}; %Hk interval data
                %call the function that generates the image
                [I_Hk] = grid_gaussian_bump(Hk, grid_values1_Hk, grid_values2_Hk, sigma,lin_interp_vals_Hk);  
                data_images{p,t,k}=I_Hk;
                end
            end   
end
end
function [ data_images ] = soft_bound_PIs( b_p_data, max_b_p_Hk, bump_func_params, res,sig)
%soft_bound_PIs generates the PIs for a set of point clouds with the
%b_p_data. Soft refers to the fact that we add three times the variance
% to the maximal values to determine our boundaries.
%   INPUTS:        -b_p_data: birth-persistence points for each of the
%                  point clouds. Each sheet corresponds to a different
%                  Betti dimension.
%                  -max_b_p_Hk: gives the maximal persistence and maximal
%                  birth time across all point clouds for each Hk. 
%                  This information is used to create the boundaries for 
%                  the persistence images.
%                  -bump_func_params: give the values for the
%                  bump function parameters which determine how fast the
%                  function leaves zero and how quickly it levels off at
%                  1. [M1,M2] where M1 controls the speed of leaving zero
%                  and M2 controls the speed of leveling off.
%                  -res: the resolution (number of pixels). we create
%                  square images with rectangular pixels.
%                  -sig: the variance of the gaussians. 
%   OUTPUT:        -data_images: the set of persistence images generated
%                  using the selected parameter values. Each sheet
%                  corresponds to the PIs generated for different Hk
%                  interval data.



[m,n,o]=size(b_p_data);
%allocate space
data_images=cell(m,n,o);
for k=1:o    
Hk_max_b=max_b_p_Hk(k,1);
Hk_max_p=max_b_p_Hk(k,2);    
%select 100 sample values of the bump function to linearly interpolate
%between to get scaling factors.
X_Hk=linspace(0,Hk_max_p,100);    
%compute the bump function values at the test points. 
[lin_interp_vals_Hk] = approx_bump(X_Hk,0,Hk_max_p,[0,0]);    
sigma=[sig,sig]; %duplicate the variance for the PIs      
%set up gridding for Hk
birth_stepsize_Hk=(Hk_max_b+3*sig)/res; %the x-width of a pixel
persistence_stepsize_Hk=(Hk_max_p+3*sig)/res; %the y-height of a pixel
grid_values1_Hk=0:birth_stepsize_Hk:(Hk_max_b+3*sig); %must be increasing from zero to max_dist
grid_values2_Hk=(Hk_max_p+3*sig):-persistence_stepsize_Hk:0; %must be decreasing from max_dist to zero
            for p=1:m
                for t=1:n
                Hk=b_p_data{p,t,k}; %Hk interval data
                [I_Hk] = grid_gaussian_bump(Hk, grid_values1_Hk, grid_values2_Hk, sigma,lin_interp_vals_Hk);  
                data_images{p,t,k}=I_Hk;
                end
            end
end
end
function [integral_image]=grid_gaussian_bump(Hk, grid_values1_Hk, grid_values2_Hk, sigma,lin_interp_vals_Hk)
%grid_gaussian_bump takes in birth-persistence data for a single point cloud, 
%pixel boundary values,the gaussian variance, and values of the bump
%function to linearly interpolate between to get the weighting values for
%each gaussian.
%Inputs:        -BPPairs is the matrix containing the (birth, persistence)
%               pairs for each interval. This comes from calling the
%               birth_persistence_coordinate funtion.
%               -grid_values1_Hk: is a vector containing the boundaries for
%               each pixel: grid_values1 is increasing and grid_values2 is
%               decreasing
%               -sigma is a 1x2 vector with the x and y standard
%               deviations.
%               -Lin_interp_vals is 2x100 matrix containing values of
%               the bump function at a range of values up to the longest
%               persistence. Will be used to identify the scaling value
%               based on using linear interpolation between the two closest
%               persistences.
%Outputs:       -Integral image: is the image computed by discreting based
%               on the values contained in grid_values and summing over all
%               thedifferent gaussians centered at each point in the
%               birth-persistence interval data.
max_bar_length=grid_values2_Hk(1);
%keyboard
[X,Y]=meshgrid(grid_values1_Hk,grid_values2_Hk);
XX=reshape(X,[],1);
YY=reshape(Y,[],1);
LIV=lin_interp_vals_Hk;
for l=1:size(Hk,1)
    b=find(Hk(l,2)<=LIV(2:end,1) & Hk(l,2)>=LIV(1:end-1,1));
    %find the closest bump function values to the persistence of the point
    %we are working with. Find an upper and lower bound.
    M=((LIV(b+1,2)-LIV(b,2))/(LIV(b+1,1)-LIV(b,1)))*(Hk(l,2)-LIV(b,1))+LIV(b,2);
    %linearly interpolate between the upper and lower persistence bound
    %weights to get the weighting value for the Gaussian.
    AA=(M)*mvncdf([XX YY], Hk(l,:), sigma);
    AA=reshape(AA,[],length(grid_values1_Hk));
    ZZ(:,:,l)=(-AA(2:end,2:end)-AA(1:end-1,1:end-1)+AA(2:end, 1:end-1)+AA(1:end-1,2:end));
end
 integral_image=sum(ZZ,3);
end


end

