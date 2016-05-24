"UPMEsfromq" <-
function(q)
{
n=ncol(q)
N=nrow(q)
s=rep(0,times=N) 
for(k in 1:N) 
if(n!=0) if(runif(1)<q[k,n]) {s[k]=1;n=n-1} 
s 
}

