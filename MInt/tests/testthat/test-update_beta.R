cat("Running tests for update_beta ")


dFile <- system.file("extdata", "x.txt", package="MInt");
rFile <- system.file("extdata", "y.txt", package="MInt");

# Test 1.
# Tests update for beta.
m <- mint(y = rFile, x = dFile, fmla = ~ feature1 + feature2)
m <- read_data(m);
m <- initialize_optimizer(m);
m <- initialize_parameters(m); 

# Clamp all parameters except beta
m$optim$clamp$P <- TRUE;
m$optim$clamp$w <- TRUE
m$optim$clamp$mu <- TRUE;

# Set all paramters except beta to their true value
P_true <- as.matrix( read.csv(system.file("extdata", "P_true.txt", package="MInt"), header=FALSE) );
w_true <- as.matrix( read.csv(system.file("extdata", "w_true.txt", package="MInt"), header=FALSE) );
b_true <- as.matrix( read.csv(system.file("extdata", "b_true.txt", package="MInt"), header=FALSE) );
# mu_true is always the zero vector.

m$param$P <- P_true;
m$param$w <- w_true;

m <- update_beta(m);

print(cor(as.vector(b_true), as.vector(m$param$beta)))
plot(as.vector(b_true), as.vector(m$param$beta));
expect_true(cor(as.vector(b_true), as.vector(m$param$beta)) > 0.95)


# End of tests.
cat("\n");