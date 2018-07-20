function integer=randi_helper(k)

% This function returns randi(k) if integer k is at least one, and
% otherwise it returns 0. The reason is that randi(0) gives an error in 
% Matlab, and this function avoids that!

if k>=1
    integer=randi(k);
else
    integer=0;
end