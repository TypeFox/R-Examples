options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test chen-x.f ...\n")

chenl.f <- function(xx) {
x <- exp(xx)
v <- xx + exp(x)
f <- (v - sqrt(v^2 + 5e-04))/2
sum (f * f)
}

chenl.g <- function(xx) {
   res <- chenl.res(xx)
   jj<-chenl.jac(xx)
   gg<- 2.0 * jj %*% as.vector(res)
   #return(gg)
}

chenl.res <- function(xx) {
   x <- exp(xx)
   v <- xx + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
}

chenl.jac <- function(xx) {
   x <- exp(xx)
   n<-length(x)
   v <- xx + exp(x)
   jj<-matrix(0.0, n, n)
   for (i in 1:n) {
     jj[i,i] <- 0.5 * (1.0 + exp(xx[i])*exp(x[i])) * (1.0-v[i]/sqrt(v[i]^2 + 5e-04))
   } 
   jj #return(jj)
}

n<-50

p0 <- log(rexp(n)) # use log to keep away from 0 and failure

system.time(ans.optx <- optimx(par=p0, fn=chenl.f, 
   control=list(all.methods=TRUE,save.failures=TRUE,maxit=2500)))[1]

print(ans.optx)

system.time(ans.optxG <- optimx(par=p0, fn=chenl.f, gr=chenl.g, 
   control=list(all.methods=TRUE,save.failures=TRUE,maxit=2500)))[1]

print(ans.optxG)

ans.optx<-transform(ans.optx, pflag="F")
ans.optxG<-transform(ans.optxG, pflag="FG")

achenl<-rbind(ans.optx, ans.optxG)

idx<-order(as.vector(achenl$method, mode="character"))
achenl<-achenl[idx, ]
achenl
 

cat("========================= end chenlog_test ======================\n")

