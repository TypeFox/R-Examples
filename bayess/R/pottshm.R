pottshm=function(ncol=2,niter=10^4,n,m=n,beta=0){
# Metropolis-Hastings sampler for a Potts model
x=matrix(sample(1:ncol,n*m,replace=TRUE),n,m)
for (i in 1:niter){
  sampl=sample(1:(n*m))
  for (k in 1:(n*m)){
    xcur=x[sampl[k]]
    a=(sampl[k]-1)%%n+1
    b=(sampl[k]-1)%/%n+1
    xtilde=sample((1:ncol)[-xcur],1)
    acpt=beta*(xneig4(x,a,b,xtilde)-xneig4(x,a,b,xcur))
    if (log(runif(1))<acpt) x[sampl[k]]=xtilde
    }}
return(x)
}
