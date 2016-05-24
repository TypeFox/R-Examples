"UPMEqfromw" <-
function(w,n)
{
N=length(w)
expa=array(0,c(N,n))
for(i in 1:N) expa[i,1]= sum(w[i:N]) 
for(i in (N-n+1):N) expa[i,N-i+1]=exp(sum(log(w[i:N]))) 
for(i in (N-2):1) 
for(z in 2:min(N-i,n)) 
{ 
expa[i,z]=w[i]*expa[i+1,z-1]+expa[i+1,z] 
} 
q=array(0,c(N,n))
for(i in N:1) q[i,1]= w[i]/expa[i,1]
for(i in N:(N-n+1)) q[i,N-i+1]=1
for(i in (N-2):1) 
for(z in 2:min(N-i,n)) 
q[i,z] = w[i]*expa[i+1,z-1]/expa[i,z]
q
}

