library(bbmle)
f <- function() {
  maxit <- 1000
  d <- data.frame(x=0:10,
                  y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
  mle2(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
   start=list(lymax=0,lhalf=0),
   data=d,
     control=list(maxit=maxit),
   parameters=list(lymax~1,lhalf~1))
}

m1 <- f()
p <- profile(m1)
## FIXME: check results (need to save in an environment-friendly way!)
print(head(as.data.frame(p)),digits=3)
