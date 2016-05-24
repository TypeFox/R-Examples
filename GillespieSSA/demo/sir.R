#===============================================================================
# Kermack-McKendric SIR model (Brown & Rothery, 1993)
#===============================================================================

# The Kermack-McKendrick SIR model is defined as
# dS/dt = -beta*N*S
# dI/dt = beta*N*S - gamma*I
# dR/dt = gamma*I
#
# This model consists of two reactions with the following per capita rates,
# transmission: beta
# recovery:     gamma

# Define parameters
parms <- c(beta=.001, gamma=.100)

# Define system
x0 <- c(S=500, I=1, R=0)                      # Initial state vector
nu <- matrix(c(-1,0,1,-1,0,1),nrow=3,byrow=T) # State-change matrix
a  <- c("beta*S*I", "gamma*I")                # Propensity vector
tf <- 100                                     # Final time
simName <- "Kermack-McKendrick SIR"

# Run the simulations
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Direct method
set.seed(2)
out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=TRUE,show.legend=TRUE) 

# Explicit tau-leap method
set.seed(2)
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 

# Binomial tau-leap method
set.seed(2)
out <- ssa(x0,a,nu,parms,tf,method="BTL",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 

# Optimized tau-leap method
set.seed(2)
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 
