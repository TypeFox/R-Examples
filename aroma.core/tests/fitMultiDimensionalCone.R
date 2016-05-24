library("aroma.core")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit genotype cone based on available methods (==packages)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
flavors <- c("sfit", "expectile")[1]
keep <- sapply(flavors, FUN=requireNamespace, quietly=TRUE)
flavors <- flavors[keep]
cat("Available fitting flavors:", paste(flavors, collapse=", "), "\n")
hasSfit <- is.element("sfit", flavors);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data (taken from the cfit.matrix() example of 'sfit')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N <- 1000

# Simulate sequences
nucleotides <- c("A", "C", "G", "T")
g <- sample(nucleotides, size=N, replace=TRUE)
ndim <- length(nucleotides)

# Simulate concentrations of allele A and allele B
X <- matrix(rexp(N), nrow=N, ncol=ndim)
colnames(X) <- nucleotides
for (nucleotide in nucleotides) {
  cc <- match(nucleotide, nucleotides)
  X[g == nucleotide, -cc] <- 0
}

# The true offset
a0 <- rep(0.3, times=ndim)

# The crosstalk matrix
A <- matrix(c(
  0.9, 0.3, 0.2, 0.1,
  0.1, 0.8, 0.1, 0.1,
  0.3, 0.4, 0.7, 0.1,
  0.1, 0.1, 0.6, 0.9
), nrow=ndim, byrow=TRUE)
A <- apply(A, MARGIN=2, FUN=function(u) u / sqrt(sum(u^2)))

# Simulate random errors on the input
xi <- matrix(rnorm(ndim*N, mean=0, sd=0.05), ncol=ndim)

# Generate the noisy crosstalk affected input data
Z <- t(a0 + A %*% t(X + xi))

# Generate the noisy observations of the latter
eps <- matrix(rnorm(ndim*N, mean=0, sd=0.05), ncol=ndim)
Y <- Z + eps

# Fit crosstalk model and calibrate data accordingly
if (hasSfit) {
  fit <- fitMultiDimensionalCone(Y, flavor="sfit")
  Yc <- backtransformMultiDimensionalCone(Y, fit=fit)
}

lim <- c(-0.5,6)
layout(matrix(c(1,2,3,0,4,5,0,0,6), nrow=3, ncol=3, byrow=TRUE))
par(mar=c(5,4,1,1)+0.1)
for (ii in 1:(ndim-1)) {
  for (jj in (ii+1):ndim) {
    cc <- c(jj,ii)
    labs <- nucleotides[cc]
    plot(Y[,cc], cex=0.8, xlim=lim, ylim=lim, xlab=labs[1], ylab=labs[2])
    if (hasSfit) {
      points(Yc[,cc], cex=0.8, col="blue")
      Mcc <- fit$M[c(1,1+cc),cc]
      class(Mcc) <- class(fit$M)
      lines(Mcc, lwd=2, col="red")
    }
  }
}
