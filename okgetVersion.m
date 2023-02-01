%%%%%
%%Function to obtain the forward echo trains
%AtmTwn: The forward echo trains£¬Mtm*NTwnn
%M0:Initial amplitude
% M number of echos 
% N number of Tw of HSIR
% J=numT2
% K=numT1
%ET1=1-sigma*exp(-Twn/T1k),length(Twn)*numT1; The first kernel N*K
%ET2=exp(-tm/T2j),length(tm)*numT2; The second kernel  M*J 
%pTlT2:the forward T1-T2 spectrum, size of numT1*numT2
    
%% main function
function AtmTwn=okgetVersion(M0,ET1,ET2,pT1T2)
[ET1_row,ET1_col]=size(ET1);
[ET2_row,ET2_col]=size(ET2);
[pT1T2_row,pT1T2_col]=size(pT1T2);
AtmTwn=zeros(ET2_row,ET1_row);
if(ET2_col~=pT1T2_row)
    disp('error£¡');
end
if (pT1T2_col~=ET1_col)
    disp('error£¡');
    return;
end
AtmTwn=M0*ET2*pT1T2*ET1';