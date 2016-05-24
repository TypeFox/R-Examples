AB.withDescalation <-function(A,B,C,D,E,DLT){
pj=DLT
P0=rep(0,ncol=length(pj))
P1=rep(0,ncol=length(pj))
Q0=rep(0,ncol=length(pj))
Q1=rep(0,ncol=length(pj))
Q2=rep(0,ncol=length(pj))
n=rep(0,ncol=length(pj))


for (j in 1:length(pj))
{
k=ifelse(C>1,c(0:(C-1)),0)
P0[j]=sum(choose(A,k)*pj[j]^k*(1-pj[j])^(A-k))
}

for (j in 1:length(pj))
{
k=ifelse(D>C,c(C:D),C)
P1[j]=sum(choose(A,k)*pj[j]^k*(1-pj[j])^(A-k))
}

for (j in 1:length(pj))
{
for (k in C:D)
{
m=ifelse((E-k)>0,c(0:(E-k)),0)
Q0[j]=sum(choose(A,k)*pj[j]^k*(1-pj[j])^(A-k)*choose(B,m)*pj[j]^m*(1-pj[j])^(B-m))
}
}

for (j in 1:length(pj))
{
for (k in 0:(C-1))
{
m=ifelse((E-k)>0,c(0:(E-k)),0)
Q1[j]=sum(choose(A,k)*pj[j]^k*(1-pj[j])^(A-k)*choose(B,m)*pj[j]^m*(1-pj[j])^(B-m))
}
}

for (j in 1:length(pj))
{
for (k in 0:(C-1))
{
m=ifelse((E-k)>(D+1-k),c((D+1-k):(E-k)),0)
Q2[j]=sum(choose(A,k)*pj[j]^k*(1-pj[j])^(A-k)*choose(B,m)*pj[j]^m*(1-pj[j])^(B-m))
}
}

Njn=(A*P0+(A+B)*Q0)/(P0+Q0)

Pstar=rep(0,ncol=length(pj))

Pstar[1]=1-P0[1]-Q0[1]

for(k in 2:length(pj))
{
temp=prod(Q2[1:(k-1)])*(1-P0[k]-Q0[k])
Pstar[1]=Pstar[1]+temp
}

Pstar[length(pj)]=prod((P0+Q0))

Pik<-function(i,k)
{
P=(P0[i]+Q0[i])*(1-P0[k]-Q0[k])*(prod((P0[1:(i-1)]+Q0[1:(i-1)])))*(prod(Q2[(i+1):(k-1)]))
return(P)
}

for (i in 2:(length(pj)-1))
{
Pstar[i]=0
for (k in (i+1):length(pj))
{
Pstar[i]=Pstar[i]+Pik(i,k)
}
}

Njik<-function(j,i,k)
{
if (j < i){N=(A*P0[j]+(A+B)*Q0[j])/(P0[j]+Q0[j])}
else if  (i<=j & j<k){N=A+B}
else if  (j==k){N=(A*(1-P0[j]-P1[j])+(A+B)*(P1[j]-Q0[j]))/(1-P0[j]-Q0[j])}
else N=0
return(N)
}

N_ij=rep(0,ncol=length(pj))
N_j=rep(0,ncol=length(pj))

Nj<-function(j)
{
for (i in 1:(length(pj)-1))
{
N_ij[i]=0
for(k in i:(length(pj)))
{
N_ij[i]=N_ij[i]+Njik(j,i,k)*Pik(i,k)
}
}
N_ij[length(pj)]=Njik(j,length(pj),length(pj))*Pstar[length(pj)]
n=sum(N_ij)
return(n)
}

for (i in 1:length(pj))
{
N_j[i]=Njn[i]*Pstar[length(pj)]+Nj(i)
}

return(N_j)
}
