# load the RNetLogo package
library(RNetLogo)

# initialization of NetLogo
nl.path <- "C:/Program Files/NetLogo 5.3/app"
NLStart(nl.path, gui=FALSE)
model.path <- "models/Sample Models/Earth Science/Fire.nlogo"
NLLoadModel(paste(nl.path,model.path,sep="/"))

# function to simulate with specific density
sim <- function(density) {
    NLCommand("set density ", density, "setup")
    NLCommand("while [any? turtles] [go]");
    ret <- NLReport("(burned-trees / initial-trees)")
    return(ret)
}

# function to perform replicated simulations
rep.sim <- function(density, rep) {
  return(
    lapply(density, function(dens) {
      replicate(rep, sim(dens))
    })
  )
}


# perform simulation with density: 1-100, stepwidth: 1, 1 replication
# and plot the result
plot(d <- seq(1,100,1),unlist(rep.sim(d,1)),
     xlab="density", ylab="percent burned")

# density: 30-60, stepwidth: 5, 10 replications
d <- seq(30,60,5)
res <- rep.sim(d,10)
boxplot(res,names=d, xlab="density", ylab="percent burned")

# density: 30-60, stepwidth: 1, 20 replications
d <- seq(30,60,1)
res <- rep.sim(d,20)
boxplot(res,names=d, xlab="density", ylab="percent burned")


# close NetLogo
NLQuit()