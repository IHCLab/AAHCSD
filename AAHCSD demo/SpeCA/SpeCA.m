function [Output]=SpeCA(I_REF,seed,SNR)
ratio=0.05;
maxx=max(max(max(I_REF))); %Chikusei Moffett DC
I_REF=I_REF/maxx;
x_3D= I_REF;
[L1 L2 M]= size(x_3D);


[Lines,Columns,L] = size(I_REF);
nl = Lines; %number of lines
nc = Columns; % number of columns
nb = L; % number of bands
HIM = x_3D;
% clear the variables that are not used.
clear x
clear Lines
clear Columns
clear Ek
clear BANDS
clear wavlen
clear w
clear L

if ratio==0.01
    ma=1;
else
    ma=5;
end

pixel=nl*nc;
nv = nl*nc;
mb = 1;
p = ma;
new_nv = round((ratio*nb-ma)*pixel*(1/(nb+mb+1)));
seednum=1;
rng(seed(seednum))
    %% compression
    [Ya Yb A B index] = Coder(HIM,ma,mb,new_nv,seednum,'gaussian' );
    
    %% add noise
    %  SNR= 25;
    
    [la1,la2,na]=size(Ya);
    Noise1 = zeros(la1,la2,na);
        sigma = sqrt(norm(Ya(:),'fro')^2/(la1*la2)/10^(SNR/10));
        temp = sigma*randn(la1,la2,na);
        Noise1 =temp;
    Ya1=Ya+Noise1;
    
    [lb1,lb2,n]=size(Yb);
    Noise2 = zeros(lb1,lb2,n);
        sigma = sqrt(norm(Yb(:),'fro')^2/(lb1*lb2)/10^(SNR/10));
        temp = sigma*randn(lb1,lb2,n);
        Noise2 =temp;
    Yb1=Yb+Noise2;
    
    [lB1,lB2,nB]=size(B);
    Noise4 = zeros(lB1,lB2,nB);
        sigma = sqrt(norm(B(:),'fro')^2/(lB1*lB2)/10^(SNR/10));
        temp = sigma*randn(lB1,lB2,nB);
        Noise4 =temp;
    B1=B+Noise4;
    
    [lindex1,lindex2,nindex]=size(index);
    Noise3 = zeros(lindex1,lindex2,nindex);  
    temp = sigma*randn(lindex1,lindex2,nindex);
    Noise3 = temp;
    index1 = index+Noise3;
    
    
    %% reconstruction
    index1 = round(index1);
    jjj=find(index1>65536);       
    index1(jjj)=65536;
    kkk=find(index1<1);
    index1(kkk)=1;
    [Xest Mest Aest] = Decoder(Ya1,Yb1,ma, mb,new_nv, A, B1, p,index1);
    Output = reshape(Xest',nl,nc,nb);
end
