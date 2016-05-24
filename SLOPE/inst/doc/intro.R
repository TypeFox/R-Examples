## ---- echo=F-------------------------------------------------------------
set.seed(12345)

## ------------------------------------------------------------------------
# Problem parameters
n = 500                         # number of observations
p = 500                         # number of variables
k = 10                          # number of variables with nonzero coefficients
amplitude = 1.2*sqrt(2*log(p))  # signal amplitude (for noise level = 1)

# Design matrix
X <- local({
  probs = runif(p, 0.1, 0.5)
  probs = t(probs) %x% matrix(1,n,2)
  X0 = matrix(rbinom(2*n*p, 1, probs), n, 2*p)
  X0 %*% (diag(p) %x% matrix(1,2,1))
})

# Coefficient and response vectors
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero)
y = X %*% beta + rnorm(n)

## ---- echo=F, results='asis'---------------------------------------------
knitr::kable(X[1:10,1:10], col.names=1:10)

## ------------------------------------------------------------------------
library(SLOPE)

result <- SLOPE(X, y)
print(result)

## ------------------------------------------------------------------------
fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
fdp(result$selected)

## ------------------------------------------------------------------------
result <- SLOPE(X, y, fdr=0.1)
fdp(result$selected)

## ------------------------------------------------------------------------
result <- SLOPE(X, y, sigma=1)
fdp(result$selected)

## ------------------------------------------------------------------------
X.orth = sqrt(n) * qr.Q(qr(X))
y.orth = X.orth %*% beta + rnorm(n)
result = SLOPE(X.orth, y.orth, lambda='bhq')
fdp(result$selected)

## ------------------------------------------------------------------------
result = SLOPE(X, y, lambda='bhq')
fdp(result$selected)

