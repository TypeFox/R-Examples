library(EMCluster, quiet = TRUE)
set.seed(1234)

x <- iris[, 1:4]
ret <- em.EM(x, nclass = 5)
plotmd(x, ret$class)

