# Test function robust_standard
library(kernDeepStackNet)
set.seed(100)
randomNumbers <- rnorm(100000)
matTest <- matrix(randomNumbers, nrow=100000, ncol=1)
matRobust <- (matTest-0) / qnorm(0.75)
resMat <- robustStandard(matTest)
stopifnot(all(abs(resMat-matRobust)<0.1))
