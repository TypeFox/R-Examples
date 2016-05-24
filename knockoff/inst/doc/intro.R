## ----, echo=F------------------------------------------------------------
set.seed(4567)

## ------------------------------------------------------------------------
# Problem parameters
n = 600          # number of observations
p = 200          # number of variables
k = 30           # number of variables with nonzero coefficients
amplitude = 3.5  # signal amplitude (for noise level = 1)

# Problem data
X = matrix(rnorm(n*p), nrow=n, ncol=p)
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero)
y.sample <- function() X %*% beta + rnorm(n)

## ------------------------------------------------------------------------
library(knockoff)

y = y.sample()
result = knockoff.filter(X, y)
print(result)

## ------------------------------------------------------------------------
fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
fdp(result$selected)

## ------------------------------------------------------------------------
result = knockoff.filter(X, y, fdr = 0.10, statistic = knockoff.stat.fs)
fdp(result$selected)

## ------------------------------------------------------------------------
my_knockoff_stat <- function(X, X_ko, y) {
  abs(t(X) %*% y) - abs(t(X_ko) %*% y)
}
result = knockoff.filter(X, y, statistic = my_knockoff_stat)
fdp(result$selected)

## ------------------------------------------------------------------------
my_lasso_stat <- function(...) knockoff.stat.lasso_signed_max(..., nlambda=10*p)
result = knockoff.filter(X, y, statistic = my_lasso_stat)
fdp(result$selected)

