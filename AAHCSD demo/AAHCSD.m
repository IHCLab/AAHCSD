%=====================================================================
% Programmer: Tzu-Hsuan Lin
% E-mail: q38091518@gs.ncku.edu.tw
% Date: 2022/02/10
% -------------------------------------------------------
% Reference:
% Chia-Hsiang Lin and Tzu-Hsuan Lin, 
% "All-addition hyperspectral compressed sensing for metasurface-driven miniaturized satellite," IEEE Transactions on Geoscience and Remote Sensing,
%  vol. 60, pp. 1-15, 2022.
%======================================================================
% A convex optimization based decoder for hyperspectral compressed sensing
% [x_3D_recons]=AAHCSD(Yh,Ym,E)
%======================================================================
%  Input
%  Ym is the upper part in (6), i.e., compressed data cube of dimension L1*L2*K.
%       (L1: vertical spatial dimension; 
%        L2: horizontal spatial dimension; 
%        K: spectral dimension.)
%  Yh is the lower part in (6), i.e., compressed data cube of dimension L1/r1*L2/r2*M.
%       (L1/r1: vertical spatial dimension; 
%        L2/r2: horizontal spatial dimension; 
%        M: spectral dimension.)
%  E is the basis matrix.
%----------------------------------------------------------------------
%  Output
%  x_3D_recons is the reconstructed hyperspectral data cube of dimension L1*L2*M.
%========================================================================

function [x_3D_recons]=AAHCSD(Ym,Yh,E)
%% parameter initialize

 lambda=1; % regularization parameter
 mu=1; % penalty parameter in ADMM
 num_iter=100; % iteration of ADMM
 [L1,L2,K]=size(Ym);
 L=L1*L2;
 [lr1,lr2,M]=size(Yh);
 [~,N]=size(E);
 s=zeros(L*N,1);
 d1=zeros(L*M,1);d2=d1;
 
 %% network setting & learning
 
degreeNET= 1; % degree of network (i.e., #(edges) per patch)
patch_size= 15;
boundary_width= 0;
[Edge,Alpha]= NetLearn(mean(Ym,3),patch_size,boundary_width,degreeNET);
maxalpha=max(Alpha); 
Alpha=Alpha/maxalpha;
NoE= length(Alpha); % number of edges
Dij=zeros(L,N,NoE);
Uij=Dij;
for ij=1:NoE,
    [idx_pi_pj{ij},PTP_plus_I_inv{ij}]= patch_idx(L1,L2,Edge(ij,:),patch_size,mu,lambda,Alpha(ij)); % for Uij update
end

%% parameter preprocessing

[Perm]=Permutation(L1/lr1,lr1,lr2);
Ym_2D=reshape(Ym,[],K)';
r=L1/lr1*L2/lr2;
quo=floor(M/K); rem=mod(M,K);
b=eye(K-rem); a=eye(rem);
a1=kron(a,ones(quo+1)); b1=kron(b,ones(quo));
a2=kron(a,ones(1,quo+1)); b2=kron(b,ones(1,quo));
D=blkdiag(a2,b2);
DtD= blkdiag(a1,b1);
DtDplusmuIm=DtD+mu*eye(M);
ILkronDtym=reshape((D'*Ym_2D),[],1);
Yh_2D=reshape(Yh,lr1*lr2,M)';
BkronIMyh=reshape((kron(Yh_2D,ones(r,1)')*inv(Perm')),[],1);

%%  ADMM

s= s_int(Yh,Ym,E,Perm); % warm start
%Steven edited

for iter=1:num_iter
    
    S=reshape(s,N,L);
    
    % update u1
    nu1=reshape((E*S),[],1)-d1;
    ILkronDtymplusnu1=ILkronDtym+mu*nu1;
    rside=reshape(ILkronDtymplusnu1,M,L);
    u1=reshape((inv(DtDplusmuIm)*rside),[],1);
    
    % update u2
    nu2=reshape((E*S),[],1)-d2;
    Im1=reshape(BkronIMyh+mu*nu2,M,L);
    Im2=Im1*Perm';
    Im2=reshape(Im2,M,r,L/r);
    Im2=permute(Im2,[1,3,2]);
    Im2=reshape(Im2,L/r*M,r);
    Im3=Im2*inv(ones(r,1)*ones(r,1)'+mu*eye(r));
    Im3=reshape(Im3,M,L/r,r);
    Im3=permute(Im3,[1 3 2]);
    Im3=reshape(Im3,M,L);
    u2=reshape(Im3*inv(Perm'),[],1);
    
    % update Uij
    Delta= repmat(S',1,1,NoE)-Dij;
    Uij= Delta; % entries not in the two patches
    
    for ij=1:NoE,
        Uij(idx_pi_pj{ij},:,ij)= PTP_plus_I_inv{ij}*((mu/lambda/Alpha(ij))*Delta(idx_pi_pj{ij},:,ij)); % entries in the two patches
    end
    
    % update s
    t1=u1+d1; t2=u2+d2;
    tij=reshape(permute(Uij+Dij,[2,1,3]),[],NoE);
    sigma=sum(tij,2);
    psitt1=reshape((E'*reshape(t1,M,L)),[],1);
    psitt2=reshape((E'*reshape(t2,M,L)),[],1);
    s=1/(2+NoE)*(psitt1+psitt2+sigma);
    
    % update dual
    S=reshape(s,N,L);
    psis=reshape((E*S),[],1);
    d1=d1+u1-psis;
    d2=d2+u2-psis;
    St=repmat(S',1,1,NoE);
    Dij=Dij+Uij-St;
    
end
x_3D_recons=reshape((E*reshape(s,N,L))',L1,L2,M);

%%
function s= s_int(Yh,Ym,E,Perm);
[a1,a2,M]= size(Yh);
[L1,L2,K]= size(Ym); L=L1*L2;
r1= L1/a1;
r2= L2/a2;
Lh= L/r1/r2;

quo=floor(M/K); rem=mod(M,K);
for i=1:rem
    band_partition(1,i)=quo+1;
end
for i=rem+1:K
    band_partition(1,i)=quo;
end

Ym_block_normalize= Ym;

v= (reshape(Ym,[],K))'*Perm';
MM= reshape(v',r1*r2,Lh,K);
MM1=reshape(sum(MM,1),Lh,K);
MM2=kron(MM1,ones(r1*r2,1));
MM2= reshape(MM2,r1*r2,Lh,K);
MM= MM./MM2;
v= reshape(MM,[],K)'*Perm;
Ym_block_normalize= reshape(v',L1,L2,K);
range_start= 1;

for i=1:K,
    range= [range_start:range_start+band_partition(i)-1];
    r=length(range);
    Yhh=reshape(Yh(:,:,range),a1,[]);
    Yhh=kron(Yhh,ones(L1/a1,L2/a2));
    Yhhh(:,:,range)=reshape(Yhh,L1,L2,r);
    Ymm(:,:,range)= repmat(Ym_block_normalize(:,:,i),1,1,r);
    range_start= range(end)+1;
end

x_3D= Yhhh.*Ymm;
s= reshape(E'*(reshape(x_3D,L1*L2,M))',[],1);

return;

function [Pi,gv,B] = Permutation(r,rows_h,cols_h)
% This function helps user automatically generate the spatial spread
% transform matrix B from the blurring kernel K, automatically permute the
% pixels in a specific order (so that our fast closed-form solutions can be
% applied), and then accordingly revise the spatial spread transform matrix.
r_square=r^2; 
Temp=ones(r,r)/(r^2);

% generate B before permutation
Hv=speye(rows_h*cols_h);
for i=1:rows_h*cols_h
    Hm=sparse(reshape(Hv(:,i),rows_h,cols_h));
    B(:,i)=sparse(reshape(kron(Hm,Temp),[],1));
end
gv=reshape(Temp,[],1);

% start permutation
Pi=sparse(zeros(r_square,r_square));
Pi=kron(ones(1,rows_h*cols_h),Pi);
Pi_temp=Pi;
vv=find(B(:,1)>0);
for j=1:length(vv)
    Pi(j,vv(j))=1;
end
for i=2:rows_h*cols_h
    vv=find(B(:,i)>0);
    Pi_new=Pi_temp;
    for j=1:length(vv)
        Pi_new(j,vv(j))=1;
    end
    Pi=[Pi;Pi_new];
end
return;

function [Edge,Alpha]= NetLearn(IM,ps,bw,degree)
% [Edge,Alpha]= NetLearn(mean(Ym,3),patch_size,boundary_width,degreeNET);
[nr,nc]= size(IM);
% find overlapping patches
temp= repmat(IM,2,2); % circulant
patches= im2col(temp(1:nr+ps-1,1:nc+ps-1),[ps,ps],'sliding'); % patches(1,:) == reshape(IM,1,[])
% find invalid patch index (at the image boundary)
all_idx= ones(nr,nc);
all_idx(1+bw:nr-bw-(ps-1),1+bw:nc-bw-(ps-1))= 0; % mark valid (interior) patches by 0
invalid_idx= find(all_idx==1);
% for each non-overlapping patch, find its most similar patch
k=0;
for i= 1+bw:ps:nr-bw-(ps-1),
    for j= 1+bw:ps:nc-bw-(ps-1),
        k= k+1;
        idx= i+(j-1)*nr;
        diff= sum(((patches(:,idx)*ones(1,nr*nc))-patches).^2,1);
        diff(idx)= inf;
        diff(invalid_idx)= inf; % do not consider boundary patches
        [dist(k,1),idx_similar]= min(diff);
        Edge(k,:)= [idx,idx_similar];
        if (degree>1),
            for edge_idx=2:degree,
                k= k+1;
                diff(idx_similar)= inf;
                [dist(k,1),idx_similar]= min(diff);
                Edge(k,:)= [idx,idx_similar];
            end
        end
    end
end
dist= sqrt(dist)+1e-8;
Alpha= dist.^(-1);
return;

function [idx_pi_pj,PTP_plus_I_inv]= patch_idx(nr,nc,edge,ps,mu,lambda,alpha)
% [idx_pi_pj{ij},PTP_plus_I_inv{ij}]= patch_idx(L1,L2,Edge(ij,:),patch_size,mu,lambda,Alpha(ij));
i= edge(1);
j= edge(2);
n= nr*nc;
yi= mod(i,nr); if (yi==0), yi=nr; end
xi= 1+((i-yi)/nr);
yj= mod(j,nr); if (yj==0), yj=nr; end
xj= 1+((j-yj)/nr);
DH = @(X) reshape(sum(im2col(X,[nr,nc],'distinct'),2),nr,nc); % sum the r^2=4 subimages in X
tempi= zeros(2*nr,2*nc);
tempi(yi:yi+ps-1,xi:xi+ps-1)= ones(ps,ps);
mapi= DH(tempi);
idx_pi= find(mapi==1);
Pi= sparse(zeros(ps*ps,n));
Pi(:,idx_pi)= speye(ps*ps);
tempj= zeros(2*nr,2*nc);
tempj(yj:yj+ps-1,xj:xj+ps-1)= ones(ps,ps);
mapj= DH(tempj);
idx_pj= find(mapj==1);
Pj= sparse(zeros(ps*ps,n));
Pj(:,idx_pj)= speye(ps*ps);
idx_pi_pj= find((or(mapi,mapj))==1);
cardI= length(idx_pi_pj);
Pi_minus_Pj= sparse(Pi-Pj);
Pi_minus_Pj_active= sparse(Pi_minus_Pj(:,idx_pi_pj));
PTP_plus_I_inv= inv(((Pi_minus_Pj_active')*Pi_minus_Pj_active)+((mu/lambda/alpha)*speye(cardI)));
return;