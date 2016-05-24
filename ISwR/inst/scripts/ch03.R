sample(1:40,5)
sample(c("H","T"), 10, replace=T)
sample(c("succ", "fail"), 10, replace=T, prob=c(0.9, 0.1))
1/prod(40:36)
prod(5:1)/prod(40:36)
1/choose(40,5)
x <- seq(-4,4,0.1)
plot(x,dnorm(x),type="l")
if (.make.epsf) dev.copy2eps(file="bellcurve.ps")
x <- 0:50
plot(x,dbinom(x,size=50,prob=.33),type="h")
if (.make.epsf) dev.copy2eps(file="binomdist.ps")
1-pnorm(160,mean=132,sd=13)
pbinom(16,size=20,prob=.5)
1-pbinom(15,size=20,prob=.5)
1-pbinom(15,20,.5)+pbinom(4,20,.5)
xbar <- 83
sigma <- 12
n <- 5
sem <- sigma/sqrt(n)
sem
xbar + sem * qnorm(0.025)
xbar + sem * qnorm(0.975)
set.seed(310367)
rnorm(10)
rnorm(10)
rnorm(10,mean=7,sd=5)
rbinom(10,size=20,prob=.5)
 ## no data sets used by exercises
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
