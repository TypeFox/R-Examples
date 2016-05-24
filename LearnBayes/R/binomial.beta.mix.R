binomial.beta.mix=function(probs,betapar,data)
{
N=length(probs)
s=data[1]; f=data[2]
post.betapar=betapar+outer(rep(1,N),data)
p=post.betapar[,1]/(post.betapar[,1]+post.betapar[,2])
m.prob=exp(dbinom(s,size=s+f,prob=p,log=TRUE)+
  dbeta(p,betapar[,1],betapar[,2],log=TRUE) -
  dbeta(p,post.betapar[,1],post.betapar[,2],log=TRUE))

post.probs=probs*m.prob/sum(probs*m.prob)
return(list(probs=post.probs,betapar=post.betapar))
}