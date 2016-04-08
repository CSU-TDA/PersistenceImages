function [ vecs ] = vecs_from_PIs( PIs)
%vecs_from_PIs reshapes the persistence images into vectors to be used in a
%machine learning task.
%   Input:      -PIs: a cell array with three dimensions containing the PIs
%                as output by the function make_PIs. Each sheet (the third
%                dimension in the array) corresponds to a different betti
%                dimension
%   Output:     -vecs: a cell array with 2 dimensions containing the
%               sets of reshaped PIs. The second dimension identifies the
%               betti dimension. 
[m,n,o]=size(PIs);

vectors=cell(m,n,o);
for i=1:m
    for j=1:n
        for k=1:o
            v=reshape(PIs{i,j,k},1,[]);
            vectors{i,j,k}=v;
        end
    end
end

vecs=cell(n,o);
for j=1:n
    for k=1:o
        vecs{j,k}=cell2mat(vectors(:,j,k));
    end
end

end

