% This is the main function to conduct the numerical simulation
%Copyright (c) 2023, Xinmin Ge, China University of Petroleum 
%All rights reserved.
clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%model for kerogen and asdorbed oil
%relaxation time for the kerogen
cT1=30;
cT2=0.5;
Sc=0.6;
%relaxation time for the adsorbed oil 
aT1=25;
aT2=1.5;
Sa=0.4;
bT1=5;
bT2=5;
dT1=5;
dT2=5;
Sb=0;
Sd=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numT1=128;%number of T1
numT2=128;%number of T2
%%%Creat the forward T1-T2 matrix
[outData fT1a fT1b fT1c fT1d fT1T fT2a fT2b fT2c fT2d fT2T T1k T2j]=newgetInversionEcho(numT1,numT2,aT1,aT2,bT1,bT2,cT1,cT2,dT1,dT2,Sa,Sb,Sc,Sd);
pT1T2=outData';% The forward T1-T2 matrix
Twn=[0.1 0.5 1 5 10 20 50 100 200 1000];%Waitting time for HSIR-ME-CPMG
N=length(Twn);
%define the sequence parameter for the ME-CPMG
te1=0.06;%echo spacing first echo train
n1=100;%number of the first echo
te2=0.2;%echo spacing first echo train
n2=200;%number of the first echo
te3=0.4;%echo spacing first echo train
n3=400;%number of the first echo
te4=0.6;%echo spacing first echo train
n4=800;%number of the first echo
te5=1.2;%echo spacing first echo train
n5=1500;%number of the first echo
%M=3000;%number of the total echos
M=n1+n2+n3+n4+n5;%total number of the echos
%tm=te:te:M*te;
tm1=te1:te1:te1*n1;
tm2=tm1:te2:te2*n2;
tm3=tm2:te3:te3*n3;
tm4=tm3:te4:te4*n4;
tm5=tm4:te5:te5*n5;
tm=[tm1 tm2 tm3 tm4 tm5];
%Creat the kernel matrics for HSIR-ME-CPMG
[T1k,T2j,ET1,ET2]=okgetT1T2matrix1(numT1,numT2,tm,Twn);
M0=100;%Total amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creat the echo trains and add the noise
AtmTwn=okgetVersion(M0,ET1,ET2,pT1T2);%obtain the pure signal
snr=100;%define the noise level
AtmTwnnoise=okaddNoise(AtmTwn,snr);%add the noise
AtmTwn=AtmTwnnoise;
%%%%%non-negative constrains for echo trains
for i=1:M
    for j=1:N
        if AtmTwn(i,j)<0
            AtmTwn(i,j)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the inversion parameters
pw=0.2;%truncated value for inversion typically 0.2
T=0.005;%threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start the 2D inversion by iterative TSVD
[v1,s1,u1,s1Inv,s2Inv,v2,s2,u2,pT1T2Inv]=T1T2inversion(ET1,AtmTwn,ET2,snr,T,pw,T1k,T2j);
[WT2,WT1]=meshgrid(T2j,T1k);
%plot the inversin result
figure color white;
pcolor(WT1,WT2,pT1T2Inv);
shading interp;
set(gca,'Xscale','log','Yscale','log','FontSize',12,'FontWeight','Demi');
xlabel('T2/ms','FontSize',12,'FontWeight','Demi');ylabel('T1/ms','FontSize',12,'FontWeight','Demi');
hold on
loglog([0.01 10000],[0.01 10000],'color',[0 0 1],'linewidth',1.0);
%plot the forward model
figure color white;
pcolor(WT1,WT2,pT1T2);
shading interp;
set(gca,'Xscale','log','Yscale','log','FontSize',12,'FontWeight','Demi');
xlabel('T2/ms','FontSize',12,'FontWeight','Demi');ylabel('T1/ms','FontSize',12,'FontWeight','Demi');
hold on
loglog([0.01 10000],[0.01 10000],'color',[0 0 1],'linewidth',1.0);
%%%%%%%%%%%%%%
%Compare the T1 and T2 spectrums
f2=zeros(numT2,1);%Inverted T2 spectrum
f1=zeros(numT1,1);%Inverted T1 spectrum
f11=zeros(numT1,1);%Forward T1 spectrum
f22=zeros(numT2,1);%Forward T2 spectrum
for i=1:numT2
    f2(i)=sum(pT1T2Inv(i,:));
end
for i=1:numT1
    f1(i)=sum(pT1T2Inv(:,i));
end
for i=1:numT2
    f22(i)=sum(pT1T2(i,:));
end
for i=1:numT1
    f11(i)=sum(pT1T2(:,i));
end

