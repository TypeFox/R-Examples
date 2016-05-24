library(Rmalschains)

claw <- function(xx) {
  x <- xx[1]
  y <- (0.46 * (dnorm(x, -1, 2/3) + dnorm(x, 1, 2/3)) +
        (1/300) * (dnorm(x, -0.5, 0.01) + dnorm(x, -1,
              0.01) + dnorm(x, -1.5, 0.01)) + (7/300) *
        (dnorm(x, 0.5, 0.07) + dnorm(x, 1, 0.07) + dnorm(x,
              1.5, 0.07)))
  return(y)
}

#use MA-LS-Chains
res.claw <- malschains(function(x) {-claw(x)}, lower=c(-3), upper=c(3), verbosity=0,
                       maxEvals=50000, control=malschains.control(popsize=50, 
                       istep=300, ls="cmaes", optimum=-5))

#use only the CMA-ES local search               
res.claw2 <- malschains(function(x) {-claw(x)}, lower=c(-3), upper=c(3), verbosity=0,
                       maxEvals=50000, control=malschains.control(ls="cmaes", 
                           lsOnly=TRUE, optimum=-5))

#use only the Simplex local search               
res.claw3 <- malschains(function(x) {-claw(x)}, lower=c(-3), upper=c(3), verbosity=0,
                       maxEvals=50000, control=malschains.control(ls="simplex", 
                           lsOnly=TRUE, optimum=-5))

res.claw
res.claw2
res.claw3

x <- seq(-3, 3,length=1000)
claw_x <- NULL
for (i in 1:length(x)) claw_x[i] <- claw(x[i])

plot(x,claw_x, type="l")
points(res.claw$sol, -res.claw$fitness, col="red")
points(res.claw2$sol, pch=3, -res.claw2$fitness, col="blue")
points(res.claw3$sol, pch=3, -res.claw3$fitness, col="green")

# run the code several times to see that the local searches run fast but do not always 
# end up in the global optimum 
