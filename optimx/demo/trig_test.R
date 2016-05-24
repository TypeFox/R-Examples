options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


##########
cat("optimx test trig.f ...\n")

trig.f <- function(x){
n <- length(x)
i <- 1:n
f <- n - sum(cos(x)) + i*(1 - cos(x)) - sin(x) 
sum(f*f)
}

p0 <- rnorm(50,sd=5)
system.time(ans.optx <- optimx(par=p0, fn=trig.f, control=list(all.methods=TRUE,save.failures=TRUE,maxit=2500)))[1]
print(ans.optx)



