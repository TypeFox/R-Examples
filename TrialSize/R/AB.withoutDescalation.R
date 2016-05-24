AB.withoutDescalation <-function(A,B,C,D,E,DLT){

pj=DLT
P0=rep(0,ncol=length(pj))
P1=rep(0,ncol=length(pj))
Q0=rep(0,ncol=length(pj))
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


Nij<-function(i,j)
{
if (j < (i+1)){N=(A*P0[j]+(A+B)*Q0[j])/(P0[j]+Q0[j])}
else if  (j == (i+1)){N=(A*(1-P0[j]-P1[j])+(A+B)*(P1[j]-Q0[j]))/(1-P0[j]-Q0[j])}
else N=0
return(N)
}


N_ij=rep(0,ncol=length(pj))
N_j=rep(0,ncol=length(pj))
Pstar=rep(0,ncol=length(pj))

Pstar[1]=1-P0[1]-Q0[1]
Pstar[length(pj)]=prod((P0+Q0))

for (i in 2:(length(pj)-1))
{
Pstar[i]=(1-P0[i+1]-Q0[i+1])*(prod(P0[1:i]+Q0[1:i]))
}



Nj<-function(j)
{
for (i in 0:(length(pj)-1))
{
N_ij[i+1]=Nij(i,j)
}
n=sum(N_ij*Pstar)
return(n)
}

for (i in 1:length(pj))
{
N_j[i]=Nj(i)
}

return(N_j)
}
