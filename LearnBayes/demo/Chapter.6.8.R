##################################################
# Section 6.8 Example of Output Analysis
##################################################
 library(LearnBayes)

 d=list(int.lo=c(-Inf,seq(66,74,by=2)),
        int.hi=c(seq(66,74,by=2), Inf),
        f=c(14,30,49,70,33,15))

 library(coda)
 library(lattice)

 start=c(70,1)
 fit=laplace(groupeddatapost,start,d)

 start=c(65,1)
 proposal=list(var=fit$var,scale=0.2)
 bayesfit=rwmetrop(groupeddatapost,proposal,start,10000,d)

 dimnames(bayesfit$par)[[2]]=c("mu","log sigma")
 xyplot(mcmc(bayesfit$par[-c(1:2000),]),col="black")

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 par(mfrow=c(2,1))
 autocorr.plot(mcmc(bayesfit$par[-c(1:2000),]),auto.layout=FALSE)
 summary(mcmc(bayesfit$par[-c(1:2000),]))
 batchSE(mcmc(bayesfit$par[-c(1:2000),]), batchSize=50)

S=readline(prompt="Type  <Return>   to continue : ")

 start=c(70,1)
 proposal=list(var=fit$var,scale=2.0)
 bayesfit=rwmetrop(groupeddatapost,proposal,start,10000,d)

 dimnames(bayesfit$par)[[2]]=c("mu","log sigma")
 sim.parameters=mcmc(bayesfit$par[-c(1:2000),])
 windows()
 xyplot(mcmc(bayesfit$par[-c(1:2000),]),col="black")

s=readline(prompt="Type  <Return>   to continue : ")

 windows()
 par(mfrow=c(2,1))
 autocorr.plot(sim.parameters,auto.layout=FALSE)
 summary(sim.parameters)
 batchSE(sim.parameters, batchSize=50)


