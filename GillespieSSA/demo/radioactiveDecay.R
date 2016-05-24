#===============================================================================
# Radioactive decay model [aka Irreversible isomerization reaction set] 
# (Gillespie, 1977)
#===============================================================================

# The radioactive decay model consists of a single species and single reaction 
# channels,
#      X --c--> 0

parms <- c(k=0.5)                        # Define parameters
x0    <- c(N=1000)                       # Initial state vector
nu    <- matrix(c(-1),nrow=1,byrow=TRUE) # State change matrix
a     <- c("k*N")                        # Propensity vector
tf    <- 20                              # Final time
simName <- "Radioactive decay model"
out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE) 
ssa.plot(out,show.title=TRUE,show.legend=FALSE)