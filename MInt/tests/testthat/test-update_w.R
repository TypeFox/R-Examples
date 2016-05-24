cat("Running tests for update_w ")

dFile <- system.file("extdata", "x.txt", package="MInt");
rFile <- system.file("extdata", "y.txt", package="MInt");

# Test 1.
# Tests update for w.
m <- mint(y = rFile, x = dFile, fmla = ~ feature1 + feature2)
m <- read_data(m);
m <- initialize_optimizer(m);
m <- initialize_parameters(m); 

# Clamp all parameters except w
m$optim$clamp$P <- TRUE;
m$optim$clamp$beta <- TRUE
m$optim$clamp$mu <- TRUE;

# Re-initialize w to random noise
m$param$w <- matrix(rnorm(nrow(m$param$w)*ncol(m$param$w)), dim(m$param$w))

# Set all paramters except w to their true value
P_true <- as.matrix( read.csv(system.file("extdata", "P_true.txt", package="MInt"), header=FALSE) );
w_true <- as.matrix( read.csv(system.file("extdata", "w_true.txt", package="MInt"), header=FALSE) );
b_true <- as.matrix( read.csv(system.file("extdata", "b_true.txt", package="MInt"), header=FALSE) );
# mu_true is always the zero vector.

m$param$P <- P_true;
m$param$beta <- b_true;

# Before the update there should be no correlation
expect_true( cor(as.vector(w_true), as.vector(m$param$w)) < 0.2 );

# Perform the update
m <- update_w(m);

# Make sure the estimate is unbiased.
me <- mean(as.vector(m$param$w) - as.vector(w_true))
std <- sd(as.vector(m$param$w) - as.vector(w_true))
expect_true(me + 2*std > 0 & me - 2*std < 0)

# Expect a high correlation between the true and estimated.
expect_true(cor(as.vector(w_true), as.vector(m$param$w)) > 0.95)


# End of tests.
cat("\n");