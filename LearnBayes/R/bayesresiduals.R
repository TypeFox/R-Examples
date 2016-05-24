bayesresiduals=function(lmfit,post,k)
{
ehat=lmfit$residuals
h=hat(model.matrix(lmfit))

prob=0*ehat
for (i in 1:length(prob))
{
z1=(k-ehat[i]/post$sigma)/sqrt(h[i])
z2=(-k-ehat[i]/post$sigma)/sqrt(h[i])
prob[i]=mean(1-pnorm(z1)+pnorm(z2))
}
return(prob)
}