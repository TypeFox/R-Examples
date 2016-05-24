cat("Running tests for intialize_beta ")

# Test 1.
# beta is initialized to the solution of an approximate log-normal fit.
# Tests the accuracy of this initialization.
#
# Read in the response and the design matrix.
dFile <- system.file("extdata", "x.txt", package="MInt");
yFile <- system.file("extdata", "y.txt", package="MInt");

m <- mint(y = yFile, x = dFile, fmla = ~ feature1 + feature2)
m <- read_data(m);
m <- initialize_optimizer(m);

b1 <- ginv(m$data$x[[1]]) %*% log(m$data$y[,1] + 1);
b2 <- ginv(m$data$x[[2]]) %*% log(m$data$y[,2] + 1);
b3 <- ginv(m$data$x[[3]]) %*% log(m$data$y[,3] + 1);

m <- initialize_beta(m);

expect_equivalent(as.matrix(m$param$beta[,1]), b1)
expect_equivalent(as.matrix(m$param$beta[,2]), b2)
expect_equivalent(as.matrix(m$param$beta[,3]), b3)

expect_false(all(as.matrix(m$param$beta[,1]) == b2))
expect_false(all(as.matrix(m$param$beta[,3]) == b2))
expect_false(all(as.matrix(m$param$beta[,1]) == b3))
expect_false(all(as.matrix(m$param$beta[,3]) == b1))
expect_false(all(as.matrix(m$param$beta[,2]) == b1))

# End of tests.
cat("\n");
