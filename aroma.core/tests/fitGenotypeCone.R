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
#set.seed(0xbeef)

N <- 1000

# Simulate genotypes
g <- sample(c("AA", "AB", "AB", "BB"), size=N, replace=TRUE)

# Simulate concentrations of allele A and allele B
X <- matrix(rexp(N), nrow=N, ncol=2)
colnames(X) <- c("A", "B")
X[g == "AA", "B"] <- 0
X[g == "BB", "A"] <- 0
X[g == "AB",] <- X[g == "AB",] / 2

# Transform noisy X
xi <- matrix(rnorm(2*N, mean=0, sd=0.05), ncol=2)
a0 <- c(0,0)+0.3
A <- matrix(c(0.9, 0.1, 0.1, 0.8), nrow=2, byrow=TRUE)
A <- apply(A, MARGIN=2, FUN=function(u) u / sqrt(sum(u^2)))
Z <- t(a0 + A %*% t(X + xi))

# Add noise to Y
eps <- matrix(rnorm(2*N, mean=0, sd=0.05), ncol=2)
Y <- Z + eps


lim <- c(-1/2,6)
xlab <- "Allele A"
ylab <- "Allele B"
plot(Y, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim)
lines(x=c(0,0,2*lim[2]), y=c(2*lim[2],0,0), col="#aaaaaa", lwd=3, lty=3)

# Different alpha sequences to illustrate the impact
alphas <- c(0.075, 0.05, 0.01, 0.03, 0.002, 0.001)

cols <- seq(from=2, to=length(alphas)+1)
legend("topright", sprintf("%.3f", alphas), col=cols, lwd=4, title="Alphas")

for (kk in seq(along=alphas)) {
  for (flavor in flavors) {
    mprintf("Flavor: %s...\n", flavor)
    fit <- fitGenotypeCone(Y, alpha=alphas[kk], flavor=flavor)
    YN <- backtransformGenotypeCone(Y, fit=fit)
    if (hasSfit) {
      sfit::radials(fit$M, lwd=3, col=cols[kk], lty=ifelse(flavor == "sfit", 1,2))
      sfit::drawApex(fit$M, col=cols[kk], pch=19, cex=2)
    }
    points(YN, col=cols[kk])
    mprintf("Flavor: %s...done\n", flavor)
  }
}
lines(x=c(0,0,2*lim[2]), y=c(2*lim[2],0,0), col="#aaaaaa", lwd=3, lty=3)
