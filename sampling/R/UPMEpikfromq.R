"UPMEpikfromq" <-function(q) 
{ 
n=ncol(q)
N=nrow(q)
pro=array(0,c(N,n))
pro[1,n]=1 
for(i in 2:N) 
for(j in 2:n) 
{ 
pro[i,j]=pro[i,j]+pro[i-1,j]*(1-q[i-1,j]) 
pro[i,j-1]=pro[i,j-1]+pro[i-1,j]*(q[i-1,j]) 
}; 
for(i in 2:N) 
{ 
pro[i,1]=pro[i,1]+pro[i-1,1]*(1-q[i-1,1]) 
} 
rowSums(pro*q) 
}

