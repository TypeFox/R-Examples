library(latentnet)

data(sampson)

badfit<-ergmm(samplike~euclidean(d=2,G=3)+rreceiver,control=ergmm.control(mle.maxit=3,burnin=0,interval=1,sample.size=1000,group.deltas=0,pilot.runs=0))

plot(badfit)

mcmc.diagnostics(badfit)
