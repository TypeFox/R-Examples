#Section 8.2.4, The coda package
library(coda)

#Section 8.3, graphical diagnoses, Example 8.1
randopar=randogibs()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 8.3, nuclear pump failure data Example 8.2
outump1=kscheck()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 8.3, heidel.diag test for Example 8.2
outump2=kscheck(heidel=T)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 8.3, noisy AR Example 8.3 
sqar(10^5)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 8.4, convergence of averages, Example 8.4
sqar(multies=TRUE)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 8.4.2, Gelman and Rubin's, Example 8.5
mump()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 8.4.3, Effective sample size
beta=randopar$beta;sigma=randopar$sigma
effectiveSize(mcmc(cbind(beta,sigma)))

S=readline(prompt="Type  <Return>   to continue : ")

#Section 8.4.4, batch means Example 8.6
beta=outump1$beta;T=length(beta)
Ts=seq(100,T,le=25)
ranj=matrix(0,ncol=2,nrow=25)
for (j in 1:25){
 aT=trunc(Ts[j]^(5/8))
 batch=rep(mean(beta),le=Ts[j]-aT)
 for (t in (aT+1):Ts[j]) batch[t-aT]=mean(beta[(t-aT):t])
 sigma=2*sqrt(sum((batch-mean(beta[1:Ts[j]]))^2)*aT/((Ts[j]-aT)*(Ts[j]-aT-1)))
 ranj[j,]=mean(beta[1:Ts[j]])+c(-sigma,+sigma)
 }

plot(Ts,ranj[,1],ylim=range(ranj),col="white",xlab="Iterations",ylab="")
polygon(c(Ts,rev(Ts)),c(ranj[,1],rev(ranj[,2])),col="grey")
lines(cumsum(beta)/(1:T),col="sienna",lwd=2)

#Section 8.5.1, Adaptive MCMC : caution, Example 8.7
adapump()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 8.5.1, Adaptive MCMC : caution, Example 8.8
sqaradap()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()
