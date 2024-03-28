function Y = insdim(X, dims)

N_dims = ndims(X);

N_insdims = length(dims);
N_newdims = N_dims + N_insdims;
insdims = (N_dims + 1):N_newdims;

% Y_dims = 1:N_dims
% 
% for i = 1:N_insdims
%     Y_dims = insel(Y_dims, insdims(i), dims(i));
% end
% 
% 
% N_dims
% N_insdims
% N_newdims 

% idx = 1:N_newdims;
% idx(dims) = [];

idx = true(1,N_newdims);
idx(dims) = false;

Y_dims = zeros(1,N_newdims);
Y_dims(idx) = 1:N_dims;
Y_dims(~idx) = insdims;

Y = permute(X,Y_dims);

end