####################################################
# Section 3.4 An Illustration of Bayesian Robustness
####################################################

 library(LearnBayes)

 quantile1=list(p=.5,x=100); quantile2=list(p=.95,x=120)
 normal.select(quantile1, quantile2)

 mu = 100
 tau = 12.16
 sigma = 15
 n = 4
 se = sigma/sqrt(4)
 ybar = c(110, 125, 140)
 tau1 = 1/sqrt(1/se^2 + 1/tau^2)
 mu1 = (ybar/se^2 + mu/tau^2) * tau1^2
 summ1=cbind(ybar, mu1, tau1)
 summ1

 tscale = 20/qt(0.95, 2)
 tscale

 par(mfrow=c(1,1))
 curve(1/tscale*dt((x-mu)/tscale,2),
   from=60, to=140, xlab="theta", ylab="Prior Density")
 curve(dnorm(x,mean=mu,sd=tau), add=TRUE, lwd=3)
 legend("topright",legend=c("t density","normal density"),
   lwd=c(1,3))

S=readline(prompt="Type  <Return>   to continue : ")

 norm.t.compute=function(ybar) {
     theta = seq(60, 180, length = 500)
     like = dnorm(theta,mean=ybar,sd=sigma/sqrt(n))
     prior = dt((theta - mu)/tscale, 2)
     post = prior * like
     post = post/sum(post)
     m = sum(theta * post)
     s = sqrt(sum(theta^2 * post) - m^2)
     c(ybar, m, s) }

summ2=t(sapply(c(110, 125, 140),norm.t.compute))
dimnames(summ2)[[2]]=c("ybar","mu1 t","tau1 t")
summ2

 cbind(summ1,summ2)

 theta=seq(60, 180, length=500)
 normpost = dnorm(theta, mu1[3], tau1)
 normpost = normpost/sum(normpost)
 windows()
 plot(theta,normpost,type="l",lwd=3,ylab="Posterior Density")
 like = dnorm(theta,mean=140,sd=sigma/sqrt(n))
 prior = dt((theta - mu)/tscale, 2)
 tpost = prior * like / sum(prior * like)  
 lines(theta,tpost)
 legend("topright",legend=c("t prior","normal prior"),lwd=c(1,3))
