##
## Test for bug when assigning NA values into an existing genotype vector
##

library(genetics)
G <- as.genotype( c("1/1","1/2") )
G[1] <- NA
