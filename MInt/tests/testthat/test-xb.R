cat("Running tests for xb ")

dFile <- system.file("extdata", "x.txt", package="MInt");
rFile <- system.file("extdata", "y.txt", package="MInt");


# Test 1.
# Tests xb calculation.
m <- mint(y = rFile, x = dFile, fmla = ~ feature1 + feature2)
m <- read_data(m);
m <- initialize_optimizer(m);
m <- initialize_parameters(m); # xb is computed in here.

xbmat <- matrix(,nrow=nrow(m$data$y),ncol=ncol(m$data$y))
xbmat[,1] <- m$data$x[[1]] %*% as.matrix(m$param$beta[,1]);
xbmat[,2] <- m$data$x[[2]] %*% as.matrix(m$param$beta[,2]);
xbmat[,3] <- m$data$x[[3]] %*% as.matrix(m$param$beta[,3]);

expect_equivalent(m$param$xb, xbmat)


# End of tests.
cat("\n");