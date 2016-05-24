#######################################################
# Section 3.5 Mixtures of Conjugate Priors
#######################################################

library(LearnBayes)

curve(.5*dbeta(x, 6, 14) + .5*dbeta(x, 14, 6), from=0, to=1,
  xlab="P", ylab="Density")

S=readline(prompt="Type  <Return>   to continue : ")

probs=c(.5,.5)
beta.par1=c(6, 14)
beta.par2=c(14, 6)
betapar=rbind(beta.par1, beta.par2)
data=c(7,3)
post=binomial.beta.mix(probs,betapar,data)
post

windows()
curve(post$probs[1]*dbeta(x,13,17)+post$probs[2]*dbeta(x,21,9),
  from=0, to=1, lwd=3, xlab="P", ylab="DENSITY")
curve(.5*dbeta(x,6,12)+.5*dbeta(x,12,6),0,1,add=TRUE)
legend("topleft",legend=c("Prior","Posterior"),lwd=c(1,3))

