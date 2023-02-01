%% 获得T1、T2布点，和T1、T2的指数矩阵
%%%%% gexinmin upc 2015/01/12
%%%% 2022/04/29 revised
%%% 2023.01.09
%%采用HSIR序列测量T1
%参数说明：

%输出参数：
    %T1k:纵向弛豫布点，个数为numT1
    %T2j:横向弛豫T2布点，个数为numT2
    %ET1=1-sigma*exp(-Twn/T1k),大小为length(Twn)*numT1; 第一个核 N*K
    %ET2=exp(-tm/T2j),大小为 length(tm)*numT2 第二个核 M*J
% 输入参数
    %numT1:纵向弛豫T1布点数
    %numT2:横向弛豫T2布点数
    %tm:一个CPMG串中回波个数 为列向量 不同回波间隔组成 0 0.2 0.4 0.6 0.8 1等
    %Twn:CPMG串数 为列向量 不同等待时间组成 1 3 5 10 20等
    %%%% 普通的IR-CPMG脉冲序列

%% 函数主体
function [T1k,T2j,ET1,ET2]=okgetT1T2matrix1(numT1,numT2,tm,Twn)
ET1=zeros(length(Twn),numT1);
ET2=zeros(length(tm),numT2);
%sigma=2;
sigma=2;
%棉籽油 0.1,10000
%对数布点
N=length(Twn);
M=length(tm);
T1k=logspace(-2,4,numT1);
T2j=logspace(-2,4,numT2);
for m=1:length(tm);
    for j=1:numT2
        ET2(m,j)=exp(-tm(m)./T2j(j));
    end
end
for n=1:length(Twn)
    for k=1:numT1
        ET1(n,k)=power(1-exp(-Twn(n)./T1k(k))*sigma,1);
        if ET1(n,k)<0;
            ET1(n,k)=0;
        end
    end
end
%%构造了两个核矩阵 K2 K1