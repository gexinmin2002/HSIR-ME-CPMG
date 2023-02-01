%%%Function used to obtain the forward T1-T2 spectrum
%%%Parameters
%%%outData:The forward 2D model, numT1*numT2
%%%imput parameters
%numT1:number for T1
%numT2:number for T2

%% Main function
function [outData fT1a fT1b fT1c fT1d fT1T fT2a fT2b fT2c fT2d fT2T T1k T2j]=newgetInversionEcho(numT1,numT2,aT1,aT2,bT1,bT2,cT1,cT2,dT1,dT2,Sa,Sb,Sc,Sd)   

outData=zeros(numT1,numT2);  %  outData=outDataa+outDatab+outDatac+outDatad
outDataa=zeros(numT1,numT2); %The first relaxation component
outDatab=zeros(numT1,numT2); %The second relaxation component
outDatac=zeros(numT1,numT2); %The third relaxation component
outDatad=zeros(numT1,numT2); %The fourth relaxation component
%%%Assuming all gaussian distributed for all relaxation components
T1k=logspace(-2,4,numT1);
T2j=logspace(-2,4,numT2);
for k=1:numT1
    if T1k(k)<aT1
        a1=k;
    end
end
if abs(aT1-T1k(a1))<abs(aT1-T1k(a1+1))
    a1=a1;
else
    a1=a1+1;
end
for j=1:numT2
    if T2j(j)<aT2
        a2=j;
    end
end
if abs(aT2-T2j(a2))<abs(aT2-T2j(a2+1))
    a2=a2;
else
    a2=a2+1;
end
da=0.1*sqrt(a1*a2);
for k=1:numT1
    for j=1:numT2
        outDataa(k,j)=(k-a1).*(k-a1)+(j-a2).*(j-a2);
    end
end
outDataa=-outDataa/2/da;
outDataa=exp(outDataa)/sqrt(2*pi)/sqrt(da);
%%eliminate the small data
for k=1:numT1
    for j=1:numT2
        if outDataa(k,j)<1e-2
            outDataa(k,j)=0;
        end
    end
end

T1k=logspace(-2,4,numT1);
T2j=logspace(-2,4,numT2);
for k=1:numT1
    if T1k(k)<bT1
        b1=k;
    end
end
if abs(bT1-T1k(b1))<abs(bT1-T1k(b1+1))
    b1=b1;
else
    b1=b1+1;
end


for j=1:numT2
    if T2j(j)<bT2
        b2=j;
    end
end
if abs(bT2-T2j(b2))<abs(bT2-T2j(b2+1))
    b2=b2;
else
    b2=b2+1;
end
db=sqrt(b1*b2);
for k=1:numT1
    for j=1:numT2
        outDatab(k,j)=20*(k-b1).*(k-b1)+(j-b2).*(j-b2);
    end
end
outDatab=-outDatab/2/db;
outDatab=exp(outDatab)/sqrt(2*pi)/sqrt(db);
for k=1:numT1
    for j=1:numT2
        if outDatab(k,j)<1e-2
            outDatab(k,j)=0;
        end
    end
end
T1k=logspace(-2,4,numT1);
T2j=logspace(-2,4,numT2);
for k=1:numT1
    if T1k(k)<cT1
        c1=k;
    end
end
if abs(cT1-T1k(c1))<abs(cT1-T1k(c1+1))
    c1=c1;
else
    c1=c1+1;
end
for j=1:numT2
    if T2j(j)<cT2
        c2=j;
    end
end
if abs(cT2-T2j(c2))<abs(cT2-T2j(c2+1))
    c2=c2;
else
    c2=c2+1;
end
dc=sqrt(c1*c2);
for k=1:numT1
    for j=1:numT2
        outDatac(k,j)=(k-c1).*(k-c1)+20*(j-c2).*(j-c2);
    end
end
 %%Notice: you should choose your parameters to control the position and the shape 
outDatac=-outDatac/2/dc;
outDatac=exp(outDatac)/sqrt(2*pi)/sqrt(dc);
for k=1:numT1
    for j=1:numT2
        if outDatac(k,j)<1e-2
            outDatac(k,j)=0;
        end
    end
end
T1k=logspace(-2,4,numT1);
T2j=logspace(-2,4,numT2);
for k=1:numT1
    if T1k(k)<dT1
        d1=k;
    end
end
if abs(dT1-T1k(d1))<abs(dT1-T1k(d1+1))
    d1=d1;
else
    d1=d1+1;
end
for j=1:numT2
    if T2j(j)<dT2
        d2=j;
    end
end
if abs(dT2-T2j(d2))<abs(dT2-T2j(d2+1))
    d2=d2;
else
    d2=d2+1;
end
dd=0.1*sqrt(d1*d2);
for k=1:numT1
    for j=1:numT2
        outDatad(k,j)=(k-d1).*(k-d1)+(j-d2).*(j-d2);
    end
end
outDatad=-outDatad/2/dd;
outDatad=exp(outDatad)/sqrt(2*pi)/sqrt(dd);
for k=1:numT1
    for j=1:numT2
        if outDatad(k,j)<1e-2
            outDatad(k,j)=0;
        end
    end
end
%%% Add the physical contrains T1 larger than T2
for k=1:numT1
    for j=1:k
        if (T1k(k)<T2j(j))&(outDataa(k,j)>0)
            outDataa(k,j)=0;
        end
    end
end

for k=1:numT1
    for j=1:k
        if (T1k(k)<T2j(j))&(outDatab(k,j)>0)
            outDatab(k,j)=0;
        end
    end
end
for k=1:numT1
    for j=1:k
        if (T1k(k)<T2j(j))&(outDatac(k,j)>0)
            outDatac(k,j)=0;
        end
    end
end
for k=1:numT1
    for j=1:k
        if (T1k(k)<T2j(j))&(outDatad(k,j)>0)
            outDatad(k,j)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%% Normalize the amplitude%%%%%%%%%%%%%%%%%%%%%
por=1;           % total amplitude
suma=0;
sumb=0;
sumc=0;
sumd=0;
fT2a=zeros(numT2,1);
fT2b=zeros(numT2,1);
fT2c=zeros(numT2,1);
fT2d=zeros(numT2,1);

fT1a=zeros(numT1,1);
fT1b=zeros(numT1,1);
fT1c=zeros(numT1,1);
fT1d=zeros(numT1,1);
fT2a=sum(outDataa)';
fT2b=sum(outDatab)';
fT2c=sum(outDatac)';
fT2d=sum(outDatad)';

fT1a=sum(outDataa');
fT1b=sum(outDatab');
fT1c=sum(outDatac');
fT1d=sum(outDatad');
suma=sum(fT2a);
sumb=sum(fT2b);
sumc=sum(fT2c);
sumd=sum(fT2d);
sumaa=sum(fT1a);
sumba=sum(fT1b);
sumca=sum(fT1c);
sumda=sum(fT1d);
for j=1:numT2;
    fT2a(j)=fT2a(j)/suma*Sa*por;
    fT2b(j)=fT2b(j)/sumb*Sb*por;
    fT2c(j)=fT2c(j)/sumc*Sc*por;
    fT2d(j)=fT2d(j)/sumd*Sd*por;
    for i=1:numT1
        outDataa(i,j)=outDataa(i,j)/suma*Sa*por;
        outDatab(i,j)=outDatab(i,j)/sumb*Sb*por;
        outDatac(i,j)=outDatac(i,j)/sumc*Sc*por;
        outDatad(i,j)=outDatad(i,j)/sumd*Sd*por;
    end
end
for i=1:numT1
    for j=1:numT2
        outData(i,j)=outDataa(i,j)+outDatab(i,j)+outDatac(i,j)+outDatad(i,j);
    end
end
for j=1:numT1
    fT1a(j)=fT1a(j)/suma*Sa*por;
    fT1b(j)=fT1b(j)/sumb*Sb*por;
    fT1c(j)=fT1c(j)/sumc*Sc*por;
    fT1d(j)=fT1d(j)/sumd*Sd*por;
end
fT2T=fT2a+fT2b+fT2c+fT2d;
fT1T=fT1a+fT1b+fT1c+fT1d;



