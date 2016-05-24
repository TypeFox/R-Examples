#Loading package
library(R0)

## In this example we simulate n=100 epidemic curves, with peak value at 150 incident cases, 
## and maximum epidemic length of 30 time units.
## Only the outbreak phase is computed. When the peak value is reached, the process is stopped 
## and another epidemic is generated.
sim.epid(epid.nb=100, GT=generation.time("gamma",c(3,1.5)), R0=1.5, 
         epid.length=30, family="poisson", peak.value=150)

# Here, a 30*100 matrix is returned. Each column is a single epidemic.
