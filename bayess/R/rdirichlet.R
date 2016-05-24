rdirichlet=function(n=1,par=rep(1,2)){
k=length(par)
mat=matrix(0,n,k)
for (i in 1:n){
sim=rgamma(k,shape=par,scale=1)
mat[i,]=sim/sum(sim)
}
mat
}
