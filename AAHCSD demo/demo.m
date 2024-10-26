clear all; close all; 
clc;
%% real hyperspectral image
addpath(genpath('SpeCA'))  
addpath(genpath('data'))  
load Chikusei;
load seed;
maxx=max(max(max(I_REF))); 
I_REF=I_REF/maxx;
x_3D= I_REF;
[L1 L2 M]= size(x_3D);
N=10;
[E]= compute_basis(x_3D,N);

%% parameters
SNR=45; % Gaussian noise: set as 25/35/45 db in paper
ratio= 4; % parameter in measurement matrix for spatial encoding
K=6; % parameter in measurement matrix for spectral encoding

%% encoding
r1=ratio; r2=ratio;
quo=floor(M/K); rem=mod(M,K);
for i=1:rem
    band_partition(1,i)=quo+1;
end
for i=rem+1:K
    band_partition(1,i)=quo;
end

range_start= 1;
for i=1:length(band_partition),
    range= [range_start:range_start+band_partition(i)-1];
    Ym(:,:,i)= sum(x_3D(:,:,range),3);
    range_start= range(end)+1;
end

for i=1:L1/r1,
    Yh1(i,:,:)= sum(x_3D(1+(i-1)*r1:i*r1,:,:),1);
end
for i=1:L2/r2,
    Yh2(:,i,:)= sum(Yh1(:,1+(i-1)*r2:i*r2,:),2);
end
Yh= Yh2;
Yh_DR= reshape((E'*(reshape(Yh,[],M))')',L1/r1,L2/r2,N);

%% add noise to simulate the transmission from the satellite to the ground station
[Yh_DR1,Ym1,E1]=addNOISE(Yh_DR,Ym,E,seed(1),SNR);

%% decoding
tic
Yh_in= reshape((E1*(reshape(Yh_DR1,[],N))')',L1/r1,L2/r2,M);
[x_3D_recons]=AAHCSD(Ym1,Yh_in,E1);
time=toc;

%% benchmark: SpeCA
[Out_SpeCA]=SpeCA(I_REF,seed,SNR);

%% Plot
band_set=[56 34 19]; % Chikusei
normColor=@(R)max(min((R-mean(R(:)))/std(R(:)),2),-2)/3+0.5;

figure;
subplot(1,3,1)
temp_show=I_REF(:,:,band_set);
temp_show=normColor(temp_show);
imshow(temp_show); 
xlabel('(a) (Ground Truth) Reference Image'); % ground-truth hyperspectral image

subplot(1,3,2)
temp_show=x_3D_recons(:,:,band_set);
temp_show=normColor(temp_show);
imshow(temp_show); 
xlabel('(b) AAHCSD Reconstructed Output Image (Proposed)'); % reconstructed image

subplot(1,3,3)
temp_show=Out_SpeCA(:,:,band_set);
temp_show=normColor(temp_show);
imshow(temp_show); 
xlabel('(C) SpeCA Reconstructed Output Image (Benchmark)'); % reconstructed image

%% Display & Quality Metrics
fprintf('AAHCSD performance:\n')
tar=x_3D_recons;
ref=I_REF;
[~,~,bands]=size(ref);
ref=reshape(ref,[],bands);
tar=reshape(tar,[],bands);
msr=mean((ref-tar).^2,1);
max2=max(tar,[],1).^2;
out.ave=mean(10*log10(max2./msr));

[rmse_total, ~, sam, ~, ~] = quality_assessment(x_3D,x_3D_recons, 0, 1/6);
sampling_rate= 100*(length(reshape(Ym,[],1))+length(reshape(Yh_DR,[],1))+length(reshape(E,[],1)))/(length(reshape(x_3D,[],1)));
fprintf('sampling_rate: %.3f (%%)\n',sampling_rate);
fprintf('PSNR: %.3f (dB)\n',out.ave);
fprintf('SAM: %.3f (¢X)\n',sam);
fprintf('RMSE: %.3f \n',rmse_total);
fprintf('TIME: %.3f (sec.)\n',time);


