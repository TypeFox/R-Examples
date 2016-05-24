# This file creates a subsets object.
# It is designed to be called by file apprex_1d.R

level1 <- c(1,2,3,4,5,6)
level2 <- c(1,2,3,4,  6)
level3 <- c(1,2,  4,  6)

"subsets.1d" <-
  list(
       level1 = level1,
       level2 = level2,
       level3 = level3
       )
names(subsets.1d[[1]]) <- paste("level.1_run",1:length(level1),sep=".")
names(subsets.1d[[2]]) <- paste("level.2_run",1:length(level2),sep=".")
names(subsets.1d[[3]]) <- paste("level.3_run",1:length(level3),sep=".")

rm(level1)
rm(level2)
rm(level3)

