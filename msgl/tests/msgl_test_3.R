library(msgl)

# warnings = errors
options(warn=2)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

## Lambda sequence
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = FALSE)

## Sparse Group lasso

# Dense x
fit1a <- msgl(x, classes, alpha = .5, lambda = lambda, standardize = FALSE)
# (Forced) Sparse x
fit1b <- msgl(x, classes, alpha = .5, lambda = lambda, sparse.data = TRUE, standardize = FALSE)

if(max(abs(fit1a$beta[[25]]-fit1b$beta[[25]])) > 1e-5) stop()

## without intercept

# Dense x
fit1a <- msgl(x, classes, alpha = .5, lambda = lambda, intercept = FALSE, standardize = FALSE)
# (Forced) Sparse x
fit1b <- msgl(x, classes, alpha = .5, lambda = lambda, intercept = FALSE, sparse.data = TRUE, standardize = FALSE)

if(max(abs(fit1a$beta[[25]]-fit1b$beta[[25]])) > 1e-5) stop()
