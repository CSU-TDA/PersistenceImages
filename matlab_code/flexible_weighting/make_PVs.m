function [ PVs ] = make_PVs(interval_data, res, sig, weight_func,params,type)
%make_PIs generates the set of persistence vectors for the PH interval data
%stored in the cell array titled interval_data. In the special case that
%all birth times for H0 features is equal to zero, the persistence image is
%essentially a 1D curve. This function utilizes only a 1D gaussian centered
%at each point instead of a 2D gaussian as is the case in the persistence
%images.
% INPUTS:   -interval_data: is a cell array containing the H0 interval data
%            for all of the point clouds in consideration.
%            Each interval set is assumed to be an nX2 matrix with the
%            first column corresponding to the birth times of features and
%            the second column is the death time of each feature. All death
%            times must be greater than the birth time and there must be a
%            finite death time for each feature. 
%            -res: the desired resolution of the vector i.e. the
%            discritization
%            -sig: is the desired variance of the Gaussians used to compute
%            the images
%            -weight_function: the name of the weighting function to be
%            called to generate the weightings for the bars. ex
%            @linear_ramp. The weight function needs to be a function only
%            of persistence. Needs to be able to accept the
%            birth_persistence coordinates as an input.
%            -params: the set of parameter values needing to be defined fpr
%            the weighting function. 
%                   -type: refers to the declaration of a soft' or 'hard' with
%            respect to the boundaries. type=1 produces hard bounds, type=2
%            produces soft boundaries on the vectors.
%OUTPUTS:   -PIs: The set of persistence vectors generated based on the
%            options specified for the provided interval data.
[ b_p_data, max_b_p_H0, problems ] = birth_persistence_coordinates(interval_data );
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
    params=[0, max(max(max_b_p_H0(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum persistence.
    type=1; %the default boundary setting is hard
elseif nargin==2
    res=res;
    sig=.5*(max(max(max_b_p_H0(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
    weight_func=@linear_ramp; %default setting is a linear weighting function
    params=[0, max(max(max_b_p_H0(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum persistence.
    type=1; %the default boundary setting is hard.
elseif nargin==1
    res=25; %default resolution is equal to 25 pixels.
    sig=.5*(max(max(max_b_p_H0(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
    weight_func=@linear_ramp; %default setting is a linear weighting function
    params=[0, max(max(max_b_p_H0(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum persistence.
    type=1; %the default boundary setting is hard.
end
    
if type==1       
    [ data_images ] = hard_bound_PVs( b_p_data, max_b_p_H0, weight_func, params, res,sig);
elseif type==2
    [ data_images ] = soft_bound_PVs( b_p_data, max_b_p_H0, weight_func, params, res,sig);
end
    PVs=data_images;
    
function [ b_p_data, max_b_p_H0, problems] = birth_persistence_coordinates(interval_data )
%birth_persitence_coordinates takes in the interval data as output by the
%duke TDA code (birth-death coordinates) and changes them into
%birth-persistence coordinates.
%INPUTS:      -interval_data: is a cell array containing the H0 interval data
%            for all of the point clouds in consideration. 
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
[m,n]=size(interval_data);
%allocate space
max_persistences=zeros(m,n);
max_birth_times=zeros(m,n);
birth_persistence=cell(m,n);
problems=[];
for i=1:n
    for j=1:m
        B=interval_data{j,i};
        %pulls the Hk interval data for the (j,i)th point cloud.
        max_persistences(j,i)=max(B(:,2)-B(:,1));
        %computes the Hk persistence (death-birth) for the (j,i)th point cloud
        max_birth_times(j,i)=max(B(:,1));
        %determines that maximal birth time for an Hk feature for the
        %(j,i)th point cloud. We will take the maximum over all of point
        %clouds to generate non-normalized Hk PIs. 
        C=B(:,2)-B(:,1);
        birth_persistence{j,i}=[B(:,1) C];
        %birth-persistence coordinates for Hk
%        end
        D=find(C<=0);
        if length(D)>0
            problems=[problems; j,i];
        elseif length(D)==0
            problems=problems;
        end

    end
end
H0_max_birth=max(max(max_birth_times));
%determine the maximum birth time of all H0 features across the point
%clouds
H0_max_persistence=max(max(max_persistences));
%determine the maximum persistence of all H0 features across the point
%clouds
max_b_p_H0=[H0_max_birth, H0_max_persistence];
b_p_data=birth_persistence;

end

function [ data_vecs ] = hard_bound_PVs( b_p_data, max_b_p_H0, weight_func, params, res,sig)
%hard_bound_PIs generates the PIs for a set of point clouds with the
%b_p_data. Hard refers to the fact that we cut the boundaries off hard at
%the maximum values. 
%   INPUTS:        -b_p_data: birth-persistence points for each of the
%                  point clouds. 
%                  -max_b_p_H0: gives the maximal persistence and maximal
%                  birth time across all point clouds for H0. 
%                  This information is used to create the boundaries for 
%                  the persistence vectors.
%                  -weight_func: the weighting function the user specified
%                  -params: the needed paramteres for the user specified
%                  weight function
%                  -res: the resolution (number of pixels).
%                  -sig: the variance of the gaussians. 
%   OUTPUT:        -data_images: the set of persistence vectors generated
%                  using the selected parameter values. 

[m,n]=size(b_p_data);
%allocate space
data_vecs=cell(m,n);
H0_max_p=max_b_p_H0(1,2);        
%set up discritization for H0
persistence_stepsize_H0=H0_max_p/res; %the y-height of a pixel
grid_values1_H0=0:persistence_stepsize_H0:H0_max_p; %must be increasing from zero to max_persistence
            for p=1:m
                for t=1:n
                H0_both=b_p_data{p,t}; %H0 interval data
                H0=H0_both(:,2); %we only need the persistence values
                %COULD CHANGE TO BE IN TERMS OF DEATH INSTEAD OF
                %PERSISTENCE, BUT GENERATING VECTORS OVER IMAGES MEANS A
                %SINGLE INPUT TO THE WEIGHTING FUNCTION
                [weights]=arrayfun(@(row) weight_func(H0_both(row,:), params), 1:size(H0_both,1))';
                %call the function that generates the vector
                [I_H0] = one_dim_gaussian_bump(H0, grid_values1_H0, sig, weights);  
                data_vecs{p,t}=I_H0;
                end
            end   
end

function [ data_vecs ] = soft_bound_PVs(  b_p_data, max_b_p_H0, weight_func, params, res,sig)
%soft_bound_PIs generates the PVs for a set of point clouds with the
%b_p_data (H0 only). Soft refers to the fact that we add three times the variance
% to the maximal values to determine our boundaries.
%   INPUTS:        -b_p_data: birth-persistence points for each of the
%                  point clouds. 
%                  -max_b_p_H0: gives the maximal persistence and maximal
%                  birth time across all point clouds for H0
%                  This information is used to create the boundaries for 
%                  the persistence vectors.
%                  -weight_func: the weighting function the user specified
%                  -params: the needed paramteres for the user specified
%                  weight function
%                  -res: the resolution (number of pixels). we create
%                  square images with rectangular pixels.
%                  -sig: the variance of the gaussians. 
%   OUTPUT:        -data_images: the set of persistence vectors generated
%                  using the selected parameter values.

[m,n]=size(b_p_data);
%allocate space
data_vecs=cell(m,n);  
H0_max_p=max_b_p_H0(1,2);    
%set up discritization for H0
persistence_stepsize_H0=(H0_max_p+3*sig)/res; %the y-height of a pixel
grid_values1_H0=0:persistence_stepsize_H0:(H0_max_p+3*sig); %must be increasing from zero to max_persistence
            for p=1:m
                for t=1:n
                H0_both=b_p_data{p,t}; %H0 interval data
                H0=H0_both(:,2); %we only need the persistence values
                %COULD CHANGE TO BE IN TERMS OF DEATH INSTEAD OF
                %PERSISTENCE, BUT GENERATING VECTORS OVER IMAGES MEANS A
                %SINGLE INPUT TO THE WEIGHTING FUNCTION
                [weights]=arrayfun(@(row) weight_func(H0_both(row,:), params), 1:size(H0_both,1))';
                %call the function that makes the vector
                [I_H0] = one_dim_gaussian_bump(H0, grid_values1_H0,sig,weights);  
                data_vecs{p,t}=I_H0;
                end
            end
end
function [integral_vec]=one_dim_gaussian_bump(BPPairs, values, sig,weights)
%grid Gaussian takes in a persistence diagram and generates the gaussian
%gridded image of that persistence diagram. We set the parameters of the
%bump function to be 1,1.
%       Inputs: -BPPairs is the matrix containing the birth, persistence
%               pairs
%               -grid_values is a vector containing the boundaries for each
%               pixel: grid_values1 is increasing
%               -sig: the variance of the one dimensional gaussians
%               -Lin_interp_vals is 2xsamples matrix containing values of
%               the bump function at a range of values up to the longest
%               bar length. Will be used to identify the scaling value
%               based on using linear approximation between the two closest
%               lengths.
%       Outputs: -Integral_vec: The vector generated by summing over 1D
%               gaussians centered at each persistence in the H0
%               information
for i=1:size(BPPairs,1)
    M=weights(i);
    AA=(M)*normcdf(values, BPPairs(i,:), sig);
    ZZ(:,i)=(AA(2:end)-AA(1:end-1));
end
 integral_vec=sum(ZZ,2);
end


end

