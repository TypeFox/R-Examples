normal.normal.mix=function(probs,normalpar,data)
{
N=length(probs)
y=data[1]; sigma2=data[2]
prior.mean=normalpar[,1]
prior.var=normalpar[,2]
post.precision=1/prior.var+1/sigma2
post.var=1/post.precision
post.mean=(y/sigma2+prior.mean/prior.var)/post.precision

m.prob=dnorm(y,prior.mean,sqrt(sigma2+prior.var))
post.probs=probs*m.prob/sum(probs*m.prob)
return(list(probs=post.probs,normalpar=cbind(post.mean,post.var)))
}