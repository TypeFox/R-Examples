
# Set an acceptable tolerance for floating point comparisons
tol <- .Machine$double.eps^0.5

# Returns TRUE for whole numbers of numeric class
expect_that(all(is.wholenumber(as.numeric(c(1, 2, 3, 4)))), is_true())

# Returns TRUE for whole numbers of integer class
expect_that(all(is.wholenumber(as.integer(c(1, 2, 3, 4)))), is_true())

# Returns FALSE on numeric inputs with decimals
expect_that(all(is.wholenumber(c(1, 2, 3, 4) + tol)), is_false())

# Throws an error on character input
expect_that(all(is.wholenumber(as.character(c(1, 2, 3, 4)))), throws_error())

