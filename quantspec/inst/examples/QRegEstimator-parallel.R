library(snowfall)

Y <- rnorm(100) # Try 2000 and parallel computation will in fact be faster.

# Compute without using snowfall capabilities
system.time(
  qRegEst1 <- qRegEstimator(Y, levels=seq(0.25,0.75,0.25), method="fn", parallel=FALSE)
)

# Set up snowfall
sfInit(parallel=TRUE, cpus=2, type="SOCK")
sfLibrary(quantreg)
sfExportAll()

# Compare how much faster the computation is when done in parallel
system.time(
  qRegEst2 <- qRegEstimator(Y, levels=seq(0.25,0.75,0.25), method="fn", parallel=TRUE)
)

sfStop()

# Compare results
V1 <- getValues(qRegEst1)
V2 <- getValues(qRegEst2)
sum(abs(V1-V2)) # Returns: [1] 0
