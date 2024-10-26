function [E]= compute_basis(x_3D,N);
X= reshape(x_3D,[],size(x_3D,3));
X= X';

[M,~] = size(X);
% d = mean(X,2);
% U = X-d*ones(1,L);
U=X;
[eV,lambda] = eig(U*U');
diag(lambda);
E = eV(:,M-N+1:end);
% Xd = C'*(X-d*ones(1,L));