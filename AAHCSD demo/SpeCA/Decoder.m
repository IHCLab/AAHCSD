% This functions implements the reconstruction step of thethe SpeCA algorithm
% corresponding to the lines 5,13 of the pseudo-code presented in the 
% algorithm 1 of the following paper: 
%
% [1] G. Martin, J. Bioucas-Dias, Hyperspectral blind reconstruction from 
%     random spectral projections, IEEE Journal of Selected Topics in Applied 
%     Earth Observations and Remote Sensing, 2016.
%
%
% ----  Input parameters ----------------------------------------------
%
% Ya    : measurements obtained with the matrix A
% Yb    : measurements obtained with the matrices B_k
% ma    : number of measurements to be performed using matrix A
% mb    : number of measurements to be performed using matrix B
% nv    : number of pixels to be measured with matrices B_k
% A     : matrix A used in the measurement process
% B     : matrices B_k used in the measuremente process
% p     : size of the subspace to be estimated
% index : indices i_k of the pixed measured by the matrices B_k 
%
%
% --------- Output parameters -------------------------------------------
% Xest : Estimated data in matrix format with the pixel signatures by columns
% Mest : Estimated matrix M that holds the estimated basis of the subspace
% S    : Estimated matrix S that holds the coefficients associated to the 
%        basis in Mest



function [Xest Mest S] = Decoder(Ya,Yb, ma,mb,nv, A,B,p,index)

np = size(Ya,2);
nb = size(B,2);

%%%%% SVD

[Mzf Mzi] = svd(Ya*Ya');
Mzf = Mzf(:,1:p);

%% Compute S
S = pinv(Mzf'*Mzf)*Mzf'*Ya;

%% Compute D
indexrep = repmat(index,mb,1);
indexrep = reshape(indexrep,1, mb*nv);

for i = 1:p
    D(:,(i-1)*nb+1:i*nb) = bsxfun(@times, B, S(i,indexrep)');
end

%% Compute M
Mest = D\Yb;

%% Recover dataset.
Mest = reshape(Mest,nb,p);
Xest = Mest*S;

end