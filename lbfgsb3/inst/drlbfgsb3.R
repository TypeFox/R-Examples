## drlbfgsb3.R
source("lbfgsb3.R")
source("grose.R")

n <- 25L
x<-rep(3,n)
u <- rep(10, n)
l <- rep(-10, n)
it <- seq(1,n,by=2)
l[it] <- 1 ## odd variables have 1 as lower bound

ans <- lbfgsb3(x, grose.f, gr=grose.g, lower = l, upper = u,
         gs=100)
print(ans)

tmp<-readline("try with numderiv")
ansn <- lbfgsb3(x, grose.f, lower = l, upper = u,
         gs=100)
print(ansn)
