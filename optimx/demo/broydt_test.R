
options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test broydt-x.f ...\n")

broydt.f <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
sum(f*f)
}

p0 <- rnorm(50, sd=1)
system.time(ans.optx <- optimx(par=p0, fn=broydt.f,control=list(maxit=25000,save.failures=TRUE,all.methods=TRUE)))[1]

print(ans.optx)



#allpar<-ans.optx$par # ans.optx is a dataframe!
#allmeth<-ans.optx$method
#nanswer<-length(allpar)

#for (i in 1:nanswer) { # not sure this makes sense
#	curmeth<-allmeth[[i]]
#	z <- sum(ans.optx$par[[i]])
#	cat(curmeth,": ")
#	print(z, digits=16)
#}
 

