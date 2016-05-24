options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

#########################################################################################
cat("optimx test sc2 ...\n")

sc2.f <- function(x){
n <- length(x)
vec <- 1:n
sum(vec * (exp(x) - x)) / 10
}

sc2.g <- function(x){
n <- length(x)
vec <- 1:n
vec * (exp(x) - 1) / 10
}

neg.sc2.f <- function(x){
n <- length(x)
vec <- 1:n
-sum(vec * (exp(x) - x)) / 10
}

neg.sc2.g <- function(x){
n <- length(x)
vec <- 1:n
-vec * (exp(x) - 1) / 10
}

p0 <- runif(500,min=-1, max=1)
system.time(ans.optxf <- optimx(par=p0, fn=sc2.f, control=list(maxit=2500,save.failures=TRUE,all.methods=TRUE)))[1]

print(ans.optxf)



system.time(neg.ans.optxf <- optimx(par=p0, fn=neg.sc2.f, 
              control=list(maxit=2500, maximize=TRUE,save.failures=TRUE, all.methods=TRUE)))[1]

print(neg.ans.optxf)


system.time(ans.optxg <- optimx(par=p0, fn=sc2.f, gr=sc2.g,
   control=list(maxit=2500,save.failures=TRUE,all.methods=TRUE)))[1]

print(ans.optxg)


system.time(neg.ans.optxg <- optimx(par=p0, fn=neg.sc2.f, gr=neg.sc2.g,
   control=list(maxit=2500, maximize=TRUE,save.failures=TRUE,all.methods=TRUE)))[1]

print(neg.ans.optxg)



