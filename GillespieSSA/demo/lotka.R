#===============================================================================
# Lotka predator-prey model (Gillespie, 1977; Kot, 2001)
#===============================================================================

# This version of the Lotka predator-prey model is given by
# dY1/dt = c1*Y1 - c2*Y1*Y2
# dY2/dt = c2*Y1*Y2 - c3*Y2
# consisting of the three reaction channels,,
#      Y1 --c1--> Y1 + Y1 
# Y1 + Y2 --c2--> Y2 + Y2 
#      Y1 --c3--> 0

# Define parameters
parms <- c(c1=10, c2=.01, c3=10)

# Define system
x0 <- c(Y1=1000, Y2=1000)                           # Initial state vector
nu <- matrix(c(+1, -1, 0, 0, 1, -1),nrow=2,byrow=T) # State-change matrix
a  <- c("c1*Y1", "c2*Y1*Y2","c3*Y2")                # Propensity vector  
tf <- 2                                             # Final time
simName <- "Lotka predator-prey model"

# Run the simulations 
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Direct method 
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="D",simName="Lotka predator-prey model",verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=TRUE,show.legend=TRUE)

# Explict tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,tau=0.002,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

# Don't run: gives wrong results
# Binomial tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="BTL",simName,f=100,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

# Optimized tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName,epsilon=0.1,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)


