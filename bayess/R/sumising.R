sumising=function(niter=10^3,numb,beta){
#simulations by Gibbs with warmin stage
S=0
x=matrix(sample(c(0,1),numb^2,replace=TRUE),ncol=numb)
for (i in 1:niter){
s=0
sampl1=sample(1:numb)
sampl2=sample(1:numb)
for (k in 1:numb){
for (l in 1:numb){
n0=xneig4(x,sampl1[k],sampl2[l],x[sampl1[k],sampl2[l]])
n1=xneig4(x,sampl1[k],sampl2[l],1-x[sampl1[k],sampl2[l]])
if (log(runif(1))<(beta*(n1-n0))){
x[sampl1[k],sampl2[l]]=1-x[sampl1[k],sampl2[l]]
n0=n1}
s=s+n0
}}
if (2*i>niter)
S=S+s
}
return(2*S/niter)
}
