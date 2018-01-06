%This code computes the Hadamard Product of two vectors or matrices
%Written by Jhelum Chakravorty

function R = hadamard_prod(vec,mat)

[a,b]=size(mat);

for i=1:a
    for j=1:b
        R(i,j) = vec(i)*mat(i,j);
    end
end





% R(1:-k+n+1,:) = epsil*mat(1:-k+n+1,:);
% R(k+n+1:2*n+1,:) = epsil*mat(k+n+1:2*n+1,:);
% R(-k+n+2:k+n,:) = mat(-k+n+2:k+n,:);