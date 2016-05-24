#===============================================================================
# Decaying-Dimerization Reaction Set (Gillespie, 2001)
#===============================================================================

# The Decaying-Dimerization Reaction Set consists of three
# species and four reaction channels,
#      S1 --c1--> 0
# S1 + S1 --c2--> S2
#      S2 --c3--> S1 + S1
#      S2 --c4--> S3

parms <- c(c1=1.0, c2=0.002, c3=0.5, c4=0.04)  # Define parameters
x0 <- c(s1=10000, s2=0, s3=0)                  # Initial state vector
nu <- matrix(c(-1, -2, +2,  0,                 # State-change matrix
                0, +1, -1, -1,
                0,  0,  0, +1),
                nrow=3,byrow=TRUE)
a  <- c("c1*s1", "c2*s1*s1", "c3*s2", "c4*s2") # Propensity vector
tf <- 10                                       # Final time
simName <- "Decaying-Dimerizing Reaction Set"

# Run simulations 
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Direct method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1) 
ssa.plot(out,show.title=TRUE,show.legend=FALSE)

# Explict tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,tau=0.003,verbose=TRUE,consoleInterval=1) 
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 

# Binomial tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="BTL",simName,verbose=TRUE,consoleInterval=1) 
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 

# Optimized tau-leap method
set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName,verbose=TRUE,consoleInterval=1) 
ssa.plot(out,show.title=FALSE,show.legend=FALSE) 
