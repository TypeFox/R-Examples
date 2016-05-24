set.seed(1001)
x <- runif(10)
y <- 1000+x+rnorm(10,sd=0.1)
d <- data.frame(x,y)

library(bbmle)
## warning
m1 = mle2(y~dnorm(a+b*x,sd=exp(s)),start=list(a=900,b=1,s=log(0.1)),
  control=list(parscale=c(1000,1,0.1)),data=d)

m2 = mle2(y~dnorm(a+b*x,sd=exp(s)),start=list(a=900,b=1,s=log(0.1)),
  control=list(parscale=c(b=1,a=1000,s=0.1)),data=d)

m3 = mle2(y~dnorm(a+b*x,sd=exp(s)),start=list(a=1,b=1,s=log(0.1)),
  method="L-BFGS-B",lower=c(a=1100,b=2,s=-Inf),data=d)

## warning
m4 = mle2(y~dnorm(a+b*x,sd=exp(s)),start=list(b=1,a=1200,s=log(0.1)),
  method="L-BFGS-B",lower=c(2,1100,0.1),data=d)

c1 = coef(m3)[c("a","b","s")]
c2 = coef(m4)[c("a","b","s")]
if (!all(abs(c1-c2)<1e-7)) stop("mismatch")
