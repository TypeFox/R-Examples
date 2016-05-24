

# Set a tolerance for numeric comparisons
tol <- .Machine$double.eps ^ 0.5


# TODO: iterate over available methods

################################################################################
# qft.interp throws an error when given unequal vectors
test_that("qft.interp throws an error when given unequal vectors", {

expect_that(qft.interp(nil = 1:9, tb = 1:10, mito = 1:10), throws_error())
expect_that(qft.interp(nil = 1:10, tb = 1:9, mito = 1:10), throws_error())
expect_that(qft.interp(nil = 1:10, tb = 1:10, mito = 1:9), throws_error())

})

################################################################################
# qft.interp throws an error when given non-numeric vectors
test_that("qft.interp throws an error when given non-numeric vectors", {

expect_that(qft.interp(letters[1:10], 1:10, 1:10), throws_error())
expect_that(qft.interp(1:10, letters[1:10], 1:10), throws_error())
expect_that(qft.interp(1:10, 1:10, letters[1:10]), throws_error())

})

################################################################################
# qft.interp warns if given any negative values, and warning mentions which
# input variable has the problem
test_that("qft.interp warns if given any negative values, and warning mentions which input variable has the problem", {

expect_that(qft.interp(nil = 0 - tol, tb = 0, mito = 0), gives_warning("nil"))
expect_that(qft.interp(nil = 0, tb = 0 - tol, mito = 0), gives_warning("tb"))
expect_that(qft.interp(nil = 0, tb = 0, mito = 0 - tol), gives_warning("mito"))

})


################################################################################
# qft.interp warns if given values > 10
test_that("qft.interp warns if given values > 10", {

expect_that(qft.interp(nil = 10 + tol, tb = 10, mito = 10), gives_warning())
expect_that(qft.interp(nil = 10, tb = 10 + tol, mito = 10), gives_warning())
expect_that(qft.interp(nil = 10, tb = 10, mito = 10 + tol), gives_warning())

})


################################################################################
# qft.interp tolerates NA values
test_that("qft.interp tolerates NA values", {
expect_that(
    is.na(qft.interp(nil = as.numeric(NA), tb = 0.75, mito = 10)), is_true()
)

expect_that(
    is.na(qft.interp(nil = 0.10, tb = as.numeric(NA), mito = 10)), is_true()
)

expect_that(
    is.na(qft.interp(nil = 0.10, tb = 0.75, mito = as.numeric(NA))), is_true()
)

})

################################################################################
test_that("qft.cens warns and censors if given values > 10", {

# qft.cens warns if given values > 10
expect_that(qft.cens(10 + tol), gives_warning())

# qft.cens censors if given values > 10
expect_that(qft.cens(10 + tol), equals(10))

})



################################################################################
# qft.interp returns a vector of acceptable results
test_that("qft.interp returns a vector of acceptable results", {

expect_that(qft.interp(1:10, 1:10, 1:10), is_a("character"))

})



################################################################################
# If nil, tb, or mito is NA, qft.interp returns an NA




