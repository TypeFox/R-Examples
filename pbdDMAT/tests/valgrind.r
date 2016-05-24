### Checks for memory leaks
library(pbdDMAT, quietly=TRUE)
init.grid()


m <- 100
n <- 100
x <- ddmatrix("rnorm", m, n)
y <- ddmatrix("rnorm", m, 1)

cov.x <- cov(x)
beta <- lm.fit(x, y)
pca <- prcomp(x)

finalize()
