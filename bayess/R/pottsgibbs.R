pottsgibbs=function(niter,numb,beta)
{
x=sample(1:4,numb^2,prob=rep(1,4),replace=TRUE)
x=matrix(x,numb,numb)
for (i in 1:niter)
{
echan1=sample(1:numb,numb,prob=rep(1,numb)/numb)
echan2=sample(1:numb,numb,prob=rep(1,numb)/numb)
for (k in 1:numb)
{
for (l in 1:numb)
{
n=rep(0,4)
n[1]=xneig4(x,echan1[k],echan2[l],1)
n[2]=xneig4(x,echan1[k],echan2[l],2)
n[3]=xneig4(x,echan1[k],echan2[l],3)
n[4]=4-n[1]-n[2]-n[3]
x[echan1[k],echan2[l]]=sample(1:4,1,prob=exp(beta*n))
}
}
}
x
}
