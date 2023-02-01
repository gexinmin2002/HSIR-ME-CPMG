%Function for the T1-T2 inversion by TSVD
%Imput parameter£º 
%snr: sinal to noise level
%iters:maximal iterations 
%ET1=1-sigma*exp(-Twn/T1k),size of length(Twn)*numT1 OR N*K The first kernel
%ET2=exp(-tm/T2j),size of length(tm)*numT2 OR M*J The second kernel
%pT1T2:The inversion result, size of numT1*numT2 OR K*J
%%Main function
function [v1,s1,u1,s1Inv,s2Inv,v2,s2,u2,pT1T2Inv]=T1T2inversion(ET1,AtmTwn,ET2,snr,T,pw,T1k,T2j)
%%%%USE SNR for the truncated position
%INITIALIZATIONS
%%%T is the inversion residual threshold 
%%%pw is the empirical parameter used in the TSVD default value of 1.5
ET1Tmp=ET1;
ET2Tmp=ET2;
AtmTwn=AtmTwn;
%%SVD for the kernel functions
[N,K]=size(ET1);
[M,J]=size(ET2);
p=ones(K,1);
pT1T2Inv=zeros(J,K);
[u1, s1, v1]=svd(ET1Tmp,'econ'); 
[u2, s2, v2]=svd(ET2Tmp,'econ'); 
[mu1,nu1]=size(u1);
[ms1,ns1]=size(s1);
[mv1,nv1]=size(v1);
[mu,nu2]=size(u2);
[ms2,ns2]=size(s2);
[mv2,nv2]=size(v2);
[led1r,led1c]=size(s1);
[led2r,led2c]=size(s2);
slInv=zeros(led1r,led1c);
s2Inv=zeros(led2r,led2c);
sum1=0;
sum2=0;
sum11=0;
sum22=0;
for i=1:led1r
        sum1=sum1+power(s1(i,i),pw);
end
for i=1:led2r
        sum2=sum2+power(s2(i,i),pw);
end
for i=1:led1r
    sum11=sum11+power(s1(i,i),pw);
    if (sum11/sum1)>=snr/(2+snr)
        index1=i;
        break
    end
end

 for i=1:led2r
    sum22=sum22+power(s2(i,i),pw);
    if (sum22/sum2)>=snr/(1+snr)
        index2=i;
        break
    end
 end
%%%%%Setting zeros for the truncated singular values
for i=index1+1:led1r
    s1(i,i)=0;
end
for i=index2+1:led2r
    s2(i,i)=0;
end
for i=1:led1r
    for j=1:led1c
        if i==j&s1(i,j)~=0
            s1Inv(i,j)=1/s1(i,j);
        end
    end
end
for i=1:led2r
    for j=1:led2c
        if i==j&s2(i,j)~=0
            s2Inv(i,j)=1/s2(i,j);
        end
    end
end
s11=s1';
for j=1:led1r
    for i=1:led1c
        if i==j&s11(i,j)~=0
            s11Inv(i,j)=1/s11(i,j);
        end
    end
end
Iter=4000;
sit=0;
pT1T2Inv=zeros(J,K);
pT1T2Inv=v2*pinv(s2)*u2'*AtmTwn*u1*pinv(s1')*v1';
newPT1T2=pT1T2Inv;
 for j=1:J
     for k=1:K
         if T1k(k)<T2j(j)
              newPT1T2(j,k)=0;
         end
     end
 end
  for j=1:J
     for k=1:K
         if pT1T2Inv(j,k)<0
             newPT1T2(j,k)=0;
         end
     end
 end
 dA=AtmTwn-(u2*s2*v2')*newPT1T2*v1*s1'*u1';%(u1*s1*v1')';
pp=min(min(newPT1T2))-1e-3;
 while pp<0
         for j=1:J
          for k=1:K
              if newPT1T2(j,k)<1e-3
                  newPT1T2(j,k)=0;
              end
          end
      end
      %%%%Physical constrains
      for j=1:J
          for k=1:K
              if T1k(k)<T2j(j)
                   newPT1T2(j,k)=0;
              end
              if T1k(k)>5000
                  newPT1T2(j,k)=0;
              end
          end
      end        
      dA=AtmTwn-(u2*s2*v2')*(newPT1T2)*(u1*s1*v1')';
       if sit<Iter
          dP=v2*pinv(s2)*u2'*dA*u1*pinv(s1')*v1';
          newPT1T2=newPT1T2+dP;
          pT1T2Inv=newPT1T2;
          pp=min(min(pT1T2Inv))-1e-3;
          sit=sit+1;
      else
          pT1T2Inv=newPT1T2;
          break;
      end
 end
 for j=1:J
          for k=1:K
              if newPT1T2(j,k)<1e-3
                  newPT1T2(j,k)=0;
              end
          end
      end
 pT1T2Inv=newPT1T2;
 pT1T2Inv=newPT1T2/sum(sum(newPT1T2));