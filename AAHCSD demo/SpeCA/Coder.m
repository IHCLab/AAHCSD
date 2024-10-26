% This function performs the measurements of the SpeCA algorithm, 
% corresponding to the lines 2,3 of the pseudo-code presented in the 
% algorithm 1 of the paper
%
% [1] G. Martin, J. Bioucas-Dias, Hyperspectral blind reconstruction from 
%     random spectral projections, IEEE Journal of Selected Topics in Applied 
%     Earth Observations and Remote Sensing, 2016.
%
%
%
% ------ Input parameters -----------------------------------------------
%
% Xc : original dataset size nl x  nc x nb
% ma : number of measurements to be performed using matrix A
% mb : number of measurements to be performed using matrix B
% nv : number of pixels to be measured with matrices B_k
% SNR: Signal-to-Noise-Ratio in dBs of the measurements
%
% type:   measurement matrix 
%
%   'gaussian'   => Random Gaussian matrix N(0,1)
%   'hadamard'   => Random Hadamard matrix {-1,1}
%   'binary'     => Random Rademacher matrix {0,1} 
%
%
%
% ------- Output parameters --------------------------------------------
%
% Ya : measurements obtained with the matrix A.
% Yb : measurements obtained with the matrices B_k.
% A  : matrix A used in the measurement process.
% B  : matrices B_k used in the measuremente process
% index: indices i_k of the pixed measured by the matrices B_k 




function [Ya Yb A B index] = Coder(Xc, ma, mb, nv, seednum,type )
% load ck.mat
% load ck1.mat
% Xc = addnoise_WITHSEED(Xc,SNR,sum(100*ck));  
%Xc = addnoise_WITHSEED(Xc,SNR,sum(100*ck));  %if ma = 1, use ck1  
load seed
rng(seed(seednum))
[nl nc nb] = size(Xc);
np = nl*nc;
Xm = reshape(Xc,np,nb)';

% s = RandStream('mt19937ar','Seed',sum(100*ck));
% s = RandStream('mt19937ar','Seed',sum(100*ck1));
% RandStream.setGlobalStream(s);
% save randstream s
if strcmp(type,'hadamard')
    A = randn(ma,nb);
    B = randn(mb*nv,nb);
    A = (A >= 0) - (A <0);
    B = (B >= 0) - (B < 0);
else
    if strcmp(type, 'binary')
         A = randn(qf,nb);
        B = randn(mb*nv,nb);
        A = (A >= 0);
        B = (B >= 0);
    else
        A = randn(ma,nb);
        B = randn(mb*nv,nb);
    end
end

Ya = A*Xm;

index = randperm(np);
index = index(1:nv);
indexrep = repmat(index,mb,1);
indexrep = reshape(indexrep,1, mb*nv);

Yb = sum(B.*Xm(:,indexrep)',2 );

end
