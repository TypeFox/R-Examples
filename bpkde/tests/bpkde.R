library(bpkde)

set.seed(0)

datasets <- data(package = "bpkde")$results[, 3]
uks <- list(dnorm, epanechnikov, rectangular)

for(d in datasets) {
  data(list = d)
  mat <- get(d)
  print(d)
  for(uk in uks) {
    f <- bpkde(mat, kernel = uk)
    print(f)
  }
}


# Quick test of M0

data(Kurtotic)
xbin <- mvlinbin(Kurtotic)

thetas <- seq(from = -pi/2, to = pi/2, length = 91)[-91]
alphas <- rbind(cos(thetas), sin(thetas))
b <- apply(alphas, 2, function(u, Y, bw) bw(drop(Y %*% u)),
           Y = Kurtotic, bw = bw.SJ)
bws <- list(alphas = alphas, lambdas = b)

cat("test M0\n")
print(M0(xbin, bpk, bandwidths = bws))


# Quick test for 3d

X <- matrix(rnorm(3000), 1000, 3)
dev.null <- mvlinbin(X)


