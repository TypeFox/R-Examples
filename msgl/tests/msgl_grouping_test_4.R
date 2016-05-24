library(msgl)

# warnings = errors
options(warn=2)

### Tests grouping

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

### Define grouping
grouping <- 1:400
grouping[1:5] <- 1
grouping[6:10] <- 2 

## Lambda sequence
lambda <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = 1, d = 25L, lambda.min = 0.1, standardize = TRUE)

## Lasso

# Test that grouping is ignored
fit1a <- msgl(x, classes, grouping = grouping, alpha = 1, lambda = lambda, standardize = TRUE)
# (Forced) Sparse x
fit1b <- msgl(x, classes, alpha = 1, lambda = lambda, standardize = TRUE)

if( sum(predict(fit1b, x)$classes != predict(fit1a, x)$classes) > 0 ) stop()
