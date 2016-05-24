library(msgl)

# warnings = errors
options(warn=2)

### Tests grouping

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

### Define grouping

grouping <- rep(1:4, 100)

## Lambda sequence
lambda <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = 0, d = 100L, lambda.min = 0.01, standardize = FALSE)
lambda1 <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = 0, d = 100L, lambda.min = 0.01, sparse.data = TRUE, standardize = FALSE)
if(max(abs(lambda-lambda1)) > 1e-5) stop()

lambda <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = .5, d = 100L, lambda.min = 0.01, standardize = FALSE)
lambda1 <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = .5, d = 100L, lambda.min = 0.01, sparse.data = TRUE, standardize = FALSE)
if(max(abs(lambda-lambda1)) > 1e-5) stop()

lambda <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = 1, d = 100L, lambda.min = 0.01, standardize = FALSE)
lambda1 <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = 1, d = 100L, lambda.min = 0.01, sparse.data = TRUE, standardize = FALSE)
if(max(abs(lambda-lambda1)) > 1e-5) stop()
