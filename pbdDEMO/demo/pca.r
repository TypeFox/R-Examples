library(pbdDEMO, quietly=TRUE)

init.grid()

# "read" in the data
n <- 1e4
p <- 250

comm.set.seed(1234, diff=T)
dx <- ddmatrix("rnorm", nrow=n, ncol=p, bldim=c(4,4), mean=100, sd=25)

# PCA with scaling, retaining only 90% of the variation
pca <- prcomp(x=dx, retx=TRUE, scale=TRUE)

prop_var <- cumsum(pca$sdev)/sum(pca$sdev)
i <- max(min(which(prop_var > 0.9)) - 1, 1)

new_dx <- pca$x[, 1:i]

print(new_dx)

comm.cat("\nNumber of columns retained:\t", i, "\n", quiet=T)
comm.cat("Percentage of columns retained:", i/dim(dx)[2], "\n\n", quiet=T)

finalize()
