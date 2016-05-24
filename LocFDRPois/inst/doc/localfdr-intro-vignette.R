## ------------------------------------------------------------------------
library("LocFDRPois")
n0 <- 900
n1 <- 100
lambda0 <- 1
lambda1 <- 10
sim_data <- c(rpois(n0, lambda0), rpois(n1, lambda1))

## ------------------------------------------------------------------------
SummarizeLocfdr(sim_data)

