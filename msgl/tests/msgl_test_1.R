library(msgl)

# warnings = errors
options(warn=2)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

## Lambda sequence
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01, standardize = FALSE)
lambda1 <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01, sparse.data = TRUE, standardize = FALSE)

if(max(abs(lambda-lambda1)) > 1e-5) stop()