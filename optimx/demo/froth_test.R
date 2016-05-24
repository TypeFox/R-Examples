options(digits=12)
if(!require("optimx"))stop("this test requires package optimx.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)


#########################################
cat("optimx test froth-x ...\n")

froth.f <- function(p){
# Freudenstein and Roth function (Broyden, Mathematics of Computation 1965, p. 577-593)
f <- rep(NA,length(p))
f[1] <- -13 + p[1] + (p[2]*(5 - p[2]) - 2) * p[2]
f[2] <- -29 + p[1] + (p[2]*(1 + p[2]) - 14) * p[2]
sum (f * f)
}

p0 <- rpois(2,10)
system.time(ans.optx <- optimx(par=p0, fn=froth.f, control=list(all.methods=TRUE,save.failures=TRUE,maxit=2500)))[1]

print(ans.optx)


#allpar<-ans.optx$par # ans.optx is a dataframe!
#allmeth<-ans.optx$method
#nanswer<-length(allpar)

#for (i in 1:nanswer) {
#	curmeth<-allmeth[[i]]
#	z <- sum(ans.optx$par[[i]])
#	cat(curmeth,": ")
#	print(z, digits=16)
#	good   <-   10.51597043896899
#	#on Windows 10.51597043896899
#	#on Linux64 10.51597043896899
#	#on Linux32 10.51597043896899
#	if(any(abs(good - z) > 1e-12)) cat("optimx test froth FAILED for method ",curmeth,"\n")
#}
 

