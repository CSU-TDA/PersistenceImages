function [ PIs ] = make_PIs(interval_data, res, sig, weight_func,params,type)
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
%            -weight_function: the name of the weighting function to be
%            called to generate the weightings for the bars. ex
%            @linear_ramp. The weight function needs to be a function only
%            of persistence. Needs to be able to accept the
%            birth_persistence coordinates as an input.
%            -params: the set of parameter values needing to be defined fpr
%            the weighting function. 
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

if nargin>6
    error('Error: too many input arguments')
elseif nargin==6
    res=res;
    sig=sig;
    weight_func=weight_func;
    params=params;
    type=type;
elseif nargin==5
    res=res;
    sig=sig;
    weight_func=weight_func;
    params=params;
    type=1;
elseif nargin==4
    error('Error: Incomplete weight function and parameter pair');
elseif nargin==3 
    res=res;
    sig=sig;
    weight_func=@linear_ramp; %default setting is a linear weighting function
    params=[0, max(max(max_b_p_Hk(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum persistence.
    type=1; %the default boundary setting is hard
elseif nargin==2
    res=res;
    sig=.5*(max(max(max_b_p_Hk(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
    weight_func=@linear_ramp; %default setting is a linear weighting function
    params=[0, max(max(max_b_p_Hk(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum persistence.
    type=1; %the default boundary setting is hard.
elseif nargin==1
    res=25; %default resolution is equal to 25 pixels.
    sig=.5*(max(max(max_b_p_Hk(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
    weight_func=@linear_ramp; %default setting is a linear weighting function
    params=[0, max(max(max_b_p_Hk(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum persistence.
    type=1; %the default boundary setting is hard.
end
    
    
if type==1       
    [ data_images ] = hard_bound_PIs( b_p_data, max_b_p_Hk, weight_func, params, res,sig);
elseif type==2
    [ data_images ] = soft_bound_PIs( b_p_data, max_b_p_Hk, weight_func, params, res,sig);
end
    PIs=data_images;
    
function [ b_p_data, max_b_p_Hk, problems] = birth_persistence_coordinates(interval_data )
%birth_persitence_coordinates takes in the interval data as output by the
%duke TDA code (birth-death coordinates) and changes them into
%birth-persistence coordinates.
%INPUTS:     -interval_data: is a cell array containing the interval data
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
        D=find(C<0);
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
function [ data_images ] = hard_bound_PIs( b_p_data, max_b_p_Hk, weight_func, params, res,sig)
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
%                  -weight_func: the weighting function the user specified
%                  -params: the needed paramteres for the user specified
%                  weight function
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
sigma=[sig,sig]; %duplicate the variance for the PIs      
%set up gridding for Hk
birth_stepsize_Hk=Hk_max_b/res; %the x-width of a pixel
persistence_stepsize_Hk=Hk_max_p/res; %the y-height of a pixel
grid_values1_Hk=0:birth_stepsize_Hk:Hk_max_b; %must be increasing from zero to max_dist
grid_values2_Hk=Hk_max_p:-persistence_stepsize_Hk:0; %must be decreasing from max_dist to zero
            for p=1:m
                for t=1:n
                Hk=b_p_data{p,t,k}; %Hk birth persistence data
                %CHANGES TO THE WIEGHT FUNCTION INPUTS HAPPEN IN THE ROW
                %BELOW
                [weights]=arrayfun(@(row) weight_func(Hk(row,:), params), 1:size(Hk,1))';
                %call the function that makes the image
                [I_Hk] = grid_gaussian_bump(Hk, grid_values1_Hk, grid_values2_Hk, sigma,weights);  
                data_images{p,t,k}=I_Hk;
                end
            end   
end
end
function [ data_images ] = soft_bound_PIs( b_p_data, max_b_p_Hk, weight_func, params, res,sig)
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
%                  -weight_func: the weighting function the user specified
%                  -params: the needed paramteres for the user specified
%                  weight function
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
sigma=[sig,sig]; %duplicate the variance for the PIs      
%set up gridding for Hk
birth_stepsize_Hk=(Hk_max_b+3*sig)/res; %the x-width of a pixel
persistence_stepsize_Hk=(Hk_max_p+3*sig)/res; %the y-height of a pixel
grid_values1_Hk=0:birth_stepsize_Hk:(Hk_max_b+3*sig); %must be increasing from zero to max_dist
grid_values2_Hk=(Hk_max_p+3*sig):-persistence_stepsize_Hk:0; %must be decreasing from max_dist to zero
            for p=1:m
                for t=1:n
                Hk=b_p_data{p,t,k}; %birth-persistence data
                %CHANGES TO THE WIEGHT FUNCTION INPUTS HAPPEN IN THE ROW
                %BELOW
                [weights]=arrayfun(@(row) weight_func(Hk(row,:), params), 1:size(Hk,1))';
                %call the funciton that makes the image
                [I_Hk] = grid_gaussian_bump(Hk, grid_values1_Hk, grid_values2_Hk, sigma,weights);  
                data_images{p,t,k}=I_Hk;
                end
            end
end
end
function [integral_image]=grid_gaussian_bump(Hk, grid_values1_Hk, grid_values2_Hk, sigma,weights)
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
%               -weights: vector containing the weighting value for each
%               interval as determined by the user specified weighting
%               function and parameters.
%Outputs:       -Integral image: is the image computed by discreting based
%               on the values contained in grid_values and summing over all
%               thedifferent gaussians centered at each point in the
%               birth-persistence interval data.
max_bar_length=grid_values2_Hk(1);
%keyboard
[X,Y]=meshgrid(grid_values1_Hk,grid_values2_Hk);
XX=reshape(X,[],1);
YY=reshape(Y,[],1);
for l=1:size(Hk,1)
    M=weights(l);
    AA=(M)*mvncdf([XX YY], Hk(l,:), sigma);
    AA=reshape(AA,[],length(grid_values1_Hk));
    ZZ(:,:,l)=(-AA(2:end,2:end)-AA(1:end-1,1:end-1)+AA(2:end, 1:end-1)+AA(1:end-1,2:end));
end
 integral_image=sum(ZZ,3);
end


end

