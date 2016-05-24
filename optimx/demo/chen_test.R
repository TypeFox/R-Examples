options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("optimx test chen-x.f ...\n")

chen.f <- function(x) {
v <- log(x) + exp(x)
f <- (v - sqrt(v^2 + 5e-04))/2
sum (f * f)
}

p0 <- rexp(50)
system.time(ans.optx <- optimx(par=p0, fn=chen.f, lower=0, control=list(all.methods=TRUE,save.failures=TRUE,maxit=2500)))[1]

print(ans.optx)

#allpar<-ans.optx$par # ans.optx is a dataframe!
#allmeth<-ans.optx$method
#nanswer<-length(allpar)

#for (i in 1:nanswer) {
#	curmeth<-allmeth[[i]]
#	z <- sum(ans.optx$par[[i]])
#	cat(curmeth,": ")
#	print(z, digits=16)
#}
 

