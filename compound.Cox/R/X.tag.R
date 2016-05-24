X.tag=function(n,p,q,s=1){
if(q+q*s>p){warning("q+q*s>p")}
X=matrix(0,n,p)
for(i in 1:n){
  for(j in 1:q){
    A_j=runif(1,min=-0.75,max=0.75)
    temp=c(j,q+s*(j-1)+1:s)
    X[i,temp]=(runif(s+1,min=-0.75,max=0.75)+A_j)/0.612 
  }
  if( (q+q*s+1) <=p ){
    X[i,(q+q*s+1):p]=(runif(p-(q+q*s),min=-1.5,max=1.5))/0.866
  }
}
X
}
