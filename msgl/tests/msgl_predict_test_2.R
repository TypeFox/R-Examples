library(msgl)

# warnings = errors
options(warn=2)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01, standardize = TRUE)

fit.qwe <- msgl(x, classes, lambda = lambda, intercept = FALSE)

res <- predict(fit.qwe, x)
if(min(colSums(res$classes != classes)) > 0) stop()

res <- predict(fit.qwe, x, sparse.data = TRUE)
if(min(colSums(res$classes != classes)) > 0) stop()
