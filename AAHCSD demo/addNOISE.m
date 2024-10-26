function [Yh_DR1,Ym1,E1]=addNOISE(Yh_DR,Ym,E,seed,SNR)

rng(seed);

[a1,a2,a3]=size(Yh_DR);
Noise = zeros(a1,a2,a3);


    sigma = sqrt(norm(Yh_DR(:),'fro')^2/(a1*a2*a3)/10^(SNR/10));
    temp = sigma*randn(a1,a2,a3);
    Noise = temp;


Yh_DR1= Yh_DR+Noise;
 
[lp1,lp2,np]=size(Ym);
Noise1 = zeros(lp1,lp2,np);

    sigma = sqrt(norm(Ym(:),'fro')^2/(lp1*lp2*np)/10^(SNR/10));
    temp = sigma*randn(lp1,lp2,np);
    Noise1 =temp;

Ym1=Ym+Noise1;

[le1,le2,ne]=size(E);
Noise2 = zeros(le1,le2,ne);

    sigma = sqrt(norm(E(:),'fro')^2/(le1*le2*ne)/10^(SNR/10));
    temp = sigma*randn(le1,le2,ne);
    Noise2 =temp;

E1=E+Noise2;