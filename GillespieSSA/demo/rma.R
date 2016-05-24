#===============================================================================
# Rosenzweig-MacArthur predator-prey model 
# (Pineda-Krch et al., 2007, Pineda-Krch, 2008)
#===============================================================================

# The Rosenzweig-MacArthur predator-prey is defined as
# dN/dt = r(1-N/K - alpha/(1+wN))NP
# dP/dt = c*alpha/(1+wN))NP
#
# This model has five reactions with the following per capita rates,
# prey birth:     b
# prey death:     d+(b-d)N/K
# predation:      alpha/(1+wN)
# predator birth: c*alpha/(1+wN)N
# predator death: g
#
# Propensity functions:
# a1 = b * N
# a2 = (d+(b-d)N/K) * N
# a3 = alpha/(1+wN) * N * P
# a4 = c*alpha/(1+wN) * N * P
# a5 = g * P

# Define parameters 
# (B in Figure 1 in Pineda-Krch et al. 2007)
parms <- c(b=2, d=1, K=1000, alpha=0.005, 
           w=0.0025, c=2, g=2)       # Parameters
x0  <- c(N=500, P=500)               # Initial state vector
nu  <- matrix(c(+1, -1, -1,  0,  0,  # State-change matrix
                 0,  0,  0, +1, -1),     
                 nrow=2,byrow=TRUE) 
a   <- c("b*N",                      # Propensity vector
         "(d+(b-d)*N/K)*N",
         "alpha/(1+w*N)*N*P",
         "c*alpha/(1+w*N)*N*P",
         "g*P")   

tf <- 10
simName <- "Rosenzweig-MacArthur predator-prey model"

# Run the simulations
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Direct method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=TRUE,show.legend=TRUE)

# Explicit tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,tau=0.01,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

# Binomial tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="BTL",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

# Optimized tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)
