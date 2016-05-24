################################
# Section 9.2.6 An Example
################################

library(LearnBayes)

 data(birdextinct)
 attach(birdextinct)
 logtime=log(time)
 plot(nesting,logtime)
 out = (logtime > 3)
 text(nesting[out], logtime[out], label=species[out], pos = 2)	

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 plot(jitter(size),logtime,xaxp=c(0,1,1))

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 plot(jitter(status),logtime,xaxp=c(0,1,1))

##### Least-squares fit

 fit=lm(logtime~nesting+size+status,data=birdextinct,x=TRUE,y=TRUE)
 summary(fit)

##### Sampling from posterior

 theta.sample=blinreg(fit$y,fit$x,5000)

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 par(mfrow=c(2,2))
 hist(theta.sample$beta[,2],main="NESTING",
  xlab=expression(beta[1]))
 hist(theta.sample$beta[,3],main="SIZE",
  xlab=expression(beta[2]))
 hist(theta.sample$beta[,4],main="STATUS",
  xlab=expression(beta[3]))
 hist(theta.sample$sigma,main="ERROR SD",
  xlab=expression(sigma))

 apply(theta.sample$beta,2,quantile,c(.05,.5,.95))

 quantile(theta.sample$sigma,c(.05,.5,.95))

S=readline(prompt="Type  <Return>   to continue : ")

###### Estimating mean extinction times

 cov1=c(1,4,0,0)
 cov2=c(1,4,1,0)
 cov3=c(1,4,0,1)
 cov4=c(1,4,1,1)
 X1=rbind(cov1,cov2,cov3,cov4)
 mean.draws=blinregexpected(X1,theta.sample)

 c.labels=c("A","B","C","D")
 windows()
 par(mfrow=c(2,2))
 for (j in 1:4)
   hist(mean.draws[,j],
      main=paste("Covariate set",c.labels[j]),xlab="log TIME")

S=readline(prompt="Type  <Return>   to continue : ")

######## Predicting extinction times

 cov1=c(1,4,0,0)
 cov2=c(1,4,1,0)
 cov3=c(1,4,0,1)
 cov4=c(1,4,1,1)
 X1=rbind(cov1,cov2,cov3,cov4)
 pred.draws=blinregpred(X1,theta.sample)

 c.labels=c("A","B","C","D")
 windows()
 par(mfrow=c(2,2))
 for (j in 1:4)
   hist(pred.draws[,j],
      main=paste("Covariate set",c.labels[j]),xlab="log TIME")

S=readline(prompt="Type  <Return>   to continue : ")

######### Model checking via posterior predictive distribution

 pred.draws=blinregpred(fit$x,theta.sample)
 pred.sum=apply(pred.draws,2,quantile,c(.05,.95))
 par(mfrow=c(1,1))
 ind=1:length(logtime)
 windows()
 matplot(rbind(ind,ind),pred.sum,type="l",lty=1,col=1,
  xlab="INDEX",ylab="log TIME")
 points(ind,logtime,pch=19)
 out=(logtime>pred.sum[2,])
 text(ind[out], logtime[out], label=species[out], pos = 4)

S=readline(prompt="Type  <Return>   to continue : ")

######### Model checking via bayes residuals

 prob.out=bayesresiduals(fit,theta.sample,2)
 windows()
 par(mfrow=c(1,1))
 plot(nesting,prob.out)
 out = (prob.out > 0.35)
 text(nesting[out], prob.out[out], label=species[out], pos = 4)	

