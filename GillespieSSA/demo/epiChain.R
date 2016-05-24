#===============================================================================
# SIRS metapopulation model (Pineda-Krch, 2008)
#===============================================================================

# The SIRS epidemiological metapopulation model is defined in Pineda-Krch (2008)

patchPopSize <- 500 # Patch size
U <- 20 # Number of patches

# Define parameters
parms <- c(beta = 0.001, # Transmission rate
          gamma = 0.1,   # Recovery rate
            rho = 0.005, # Loss of immunity rate
        epsilon = 0.01,  # Proportion inter-patch transmissions
              N = patchPopSize) # Patch population size (constant)

# Create the named initial state vector for the U-patch system. The structure of
# x0 is as follows (assuming a patchsize of 500 individuals),
# x0 <- c(499,1,500,0,500,0 ... )
x0 <- c(S1=(patchPopSize-1), I1=1)
if (U>1) {
  for (i in (seq(2,U))) {
    patch <- c(patchPopSize, 0)
    names(patch) <- c(paste("S",i,sep=""), paste("I",i,sep="")) 
    x0 <- c(x0, patch) 
  }
}

# State-change matrix
# Define the state change matix for a single patch
#      Reaction 1   2   3   4     State
nu <- matrix(c(-1, -1,  0, +1,  # S
               +1, +1, -1,  0), # I
             nrow=2,byrow=TRUE)

# Create the propensity vector
a <- NULL
for (patch in (seq(U))) {
  i <- patch            # Intra-patch index
  if (patch==1) j <- U  # Inter-patch index    
  else j <- patch-1

  # Construct the propensity functions for the current patch   # Reaction:
  a_patch <- c(paste("(1-epsilon)*beta*S",i,"*I",i,"",sep=""), # 1. Intra-patch infection
               paste("epsilon*beta*S",i,"*I",j,"",sep=""),     # 2. Inter-patch infection
               paste("gamma*I",i,"",sep=""),                   # 3. Recovery from  infection
               paste("rho*(N-S",i,"-I",i,")",sep=""))          # 4. Loss of immunity
  a <- c(a, a_patch)
} # for()

# Run the simulations
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Run one realization using the Direct method
set.seed(1)
out <- ssa(x0, a, nu, parms, tf=50, method = "D", simName="SIRS metapopulation model", verbose=TRUE, consoleInterval=1)
ssa.plot(out,by=5,show.title=TRUE,show.legend=FALSE)

# Run one realization using the Explict tau-leap  method
set.seed(1)
out <- ssa(x0, a, nu, parms, tf=50, method = "ETL", simName="SIRS metapopulation model", verbose=TRUE, consoleInterval=1)
ssa.plot(out,by=5,show.title=FALSE,show.legend=FALSE)

# Run one realization using the Binomial tau-leap method
set.seed(1)
out <- ssa(x0, a, nu, parms, tf=50, method = "BTL", simName="SIRS metapopulation model", verbose=TRUE, consoleInterval=1)
ssa.plot(out,by=5,show.title=FALSE,show.legend=FALSE)

# Run one realization using the Optimized tau-leap method
set.seed(1)
out <- ssa(x0, a, nu, parms, tf=50, method = "OTL", simName="SIRS metapopulation model", hor=rep(2,length(x0)), verbose=TRUE, consoleInterval=1)
ssa.plot(out,by=5,show.title=FALSE,show.legend=FALSE)
