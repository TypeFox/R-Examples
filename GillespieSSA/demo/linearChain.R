#===============================================================================
# Linear Chain System (Cao et al., 2004)
#===============================================================================

# The Linear Chain System consists of M chain reactions with
# M+1 species as follows:
# S_1 --c1--> S_2 --c2-->...--cM--> S_(M+1)

# Rate parameter
parms <- c(c=1)

# Number of chain reactions
M <- 50

# Initial state vector
x0 <- c(1000, rep(0,M)) 
names(x0) <- paste("x",seq(M+1),sep="") 
  
# State-change matrix
nu <- matrix(rep(0,(M*(M+1))),ncol=M)
diag(nu) <- -1
diag(nu[2:M,]) <- +1
nu[M+1,M] <- +1

# Propensity vector
a <- as.vector(paste("c*x",seq(M),"",sep=""))

tf <- 5 # Final time
simName <- "Linear Chain System"

# Run the simulations
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Direct method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=TRUE,show.legend=FALSE)

# Explict tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,tau=0.1,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

# Binomial tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="BTL",simName,f=50,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

# Optimal tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName,verbose=TRUE,consoleInterval=1)
ssa.plot(out,show.title=FALSE,show.legend=FALSE)
