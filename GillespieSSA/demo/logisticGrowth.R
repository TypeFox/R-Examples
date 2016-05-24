#===============================================================================
# Logistic growth (Pearl-Verhulst model) (Kot, 2001)
#===============================================================================

# The logistic growth model is given by dN/dt=rN(1-N/K) 
# where N is the number (density) of indviduals at time t, K 
# is the carrying capacity of the population, r is the 
# intrinsic growth rate of the population. We assume r=b-d 
# where b is the per capita p.c. birth rate and d is the 
# p.c. death rate. 
#
# This model consists of two reaction channels,
# N ---b--->  N + N
# N ---d'---> 0
# where d'=d+(b-d)N/K. The propensity functions are a_1=bN 
# and a_2=d'N.

parms <- c(b=2, d=1, K=1000)      # Parameters
x0 <- c(N=500)                    # Initial state vector
nu <- matrix(c(+1,-1),ncol=2)     # State-change matrix
a  <- c("b*N", "(d+(b-d)*N/K)*N") # Propensity vector
tf <- 10                          # Final time
simName <- "Logistic growth" 

# Run the simulations
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Direct method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=TRUE,show.legend=FALSE) 
 
# Explict tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,tau=0.03,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 

# Run the simulation using the Binomial tau-leap method 
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="BTL",simName,f=5,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 

# Optimized tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)
