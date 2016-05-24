## set up a data frame for prediction

set.seed(1001)
f = factor(rep(letters[1:4],each=20))
x = runif(80)
u = rnorm(4)
y = rnorm(80,mean=2+x*(3+u[f]),sd=0.1)
dat = data.frame(f,x,y)

## fit a model ... could easily do by lm() but want to
##   demonstrate the problem

library(bbmle)
m1 = mle2(y~dnorm(a+b*x,sd=exp(logs)),parameters=list(b~f),data=dat,
  start=list(a=0,b=2,logs=-3))

## data frame for prediction
pp0 = expand.grid(x=seq(0,1,length=11),
  f=levels(dat$f))

## combine frame and model data: have to keep the model data
##  around, because it contain other information needed for
##  prediction.

nrow(predict(m1,pp0))



