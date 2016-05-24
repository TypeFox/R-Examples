"writesample" <-
function(n,N) 
{
if(n==N) samples=rep(1,times=N) 
else{ 
x=numeric(N) 
row=1 
for(i in (n+1):N) row=row*i 
k=1 
for(i in 1:(N-n)) k=k*i 
row=row/k 
samples=matrix(0,row,N) 
k=1 
sol=0 
x[1]=-1 
while(k<=N && k>0) 
{ 
while(x[k]<1) 
{ 
x[k]=x[k]+1 
s=0 
for(i in 1:N) 
s=s+x[i] 
if(s==n) { 
sol=sol+1 
samples[sol,]=x 
} 
else 
if(k<N){k=k+1 
x[k]=-1 
} 
} 
k=k-1 
} 
} 
samples
}

