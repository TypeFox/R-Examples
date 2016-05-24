cat("Running tests for read_data ")

dFile <- system.file("extdata", "x.txt", package="MInt");
rFile <- system.file("extdata", "y.txt", package="MInt");

x <- read.table(file=dFile, header=TRUE, row.names="Observations")
y <- read.table(file=rFile, header=TRUE, row.names="Observations")

# Test 1.
# Tests to see if read_data sucessfully reads response specific design matrices.
m <- mint(y = rFile, x = dFile, fmla = ~ feature1 + feature2)
m <- read_data(m);

expect_true(all(m$data$xd[[1]] == x))
expect_true(all(m$data$xd[[2]] == x))
expect_true(all(m$data$xd[[3]] == x))

# Test 2.
# Tests to see if read_data sucessfully reads in the response matrix.
expect_true(all(m$data$y == y))

# Test 3.
# Tests to see if read_data sucessfully reads in only the variables specified by
# the formula
m <- mint(y = rFile, x = dFile, fmla = ~ feature1 )
m <- read_data(m);

expect_true(all(m$data$xd[[1]] == as.data.frame(x[,1])))
expect_true(all(m$data$xd[[2]] == as.data.frame(x[,1])))
expect_true(all(m$data$xd[[3]] == as.data.frame(x[,1])))

# End of tests.
cat("\n");
