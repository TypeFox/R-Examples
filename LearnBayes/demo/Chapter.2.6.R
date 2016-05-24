########################
# Section 2.6 Prediction
########################

 library(LearnBayes)

 p=seq(0.05, 0.95, by=.1)
 prior = c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
 prior=prior/sum(prior)
 m=20; ys=0:20
 pred=pdiscp(p, prior, m, ys)
 cbind(0:20,pred)

 ab=c(3.26, 7.19)
 m=20; ys=0:20
 pred=pbetap(ab, m, ys)

 p=rbeta(1000,3.26, 7.19)

 y = rbinom(1000, 20, p)

 table(y)

 freq=table(y)
 ys=as.integer(names(freq))
 predprob=freq/sum(freq)
 plot(ys,predprob,type="h",xlab="y",
   ylab="Predictive Probability")

 dist=cbind(ys,predprob)

 covprob=.9
 discint(dist,covprob)
