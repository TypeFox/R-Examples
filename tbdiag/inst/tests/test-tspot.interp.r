


# Set a tolerance for numeric comparisons
tol <- .Machine$double.eps ^ 0.5


# TODO: iterate over available methods

################################################################################
# tspot.interp throws an error when given unequal vectors
test_that("tspot.interp throws an error when given unequal vectors", {

expect_that(tspot.interp(nil = 1:9, 
                         panel.a = 1:10, 
                         panel.b = 1:10, 
                         mito = 1:10), 
            throws_error()
)

expect_that(tspot.interp(nil = 1:10, 
                         panel.a = 1:9, 
                         panel.b = 1:10, 
                         mito = 1:10), 
            throws_error()
)

expect_that(tspot.interp(nil = 1:10, 
                         panel.a = 1:10, 
                         panel.b = 1:9, 
                         mito = 1:10), 
            throws_error()
)

expect_that(tspot.interp(nil = 1:10, 
                         panel.a = 1:10, 
                         panel.b = 1:10, 
                         mito = 1:9), 
            throws_error()
)

})

################################################################################
# tspot.interp throws an error when given non-numeric vectors
test_that("tspot.interp throws an error when given non-numeric vectors", {

expect_that(tspot.interp(letters[1:10], 1:10, 1:10), throws_error())
expect_that(tspot.interp(1:10, letters[1:10], 1:10), throws_error())
expect_that(tspot.interp(1:10, 1:10, letters[1:10]), throws_error())

})

################################################################################
# tspot.interp warns if given any negative values, and warning mentions which
# input variable has the problem
test_that("tspot.interp warns if given any negative values, and warning mentions which input variable has the problem", {

expect_that(tspot.interp(nil = 0 - tol, 
                         panel.a = 0, 
                         panel.b = 0, 
                         mito = 0), 
            gives_warning("nil")
)

expect_that(tspot.interp(nil = 0, 
                         panel.a = 0 - tol, 
                         panel.b = 0, 
                         mito = 0), 
            gives_warning("panel.a")
)

expect_that(tspot.interp(nil = 0, 
                         panel.a = 0, 
                         panel.b = 0 - tol, 
                         mito = 0), 
            gives_warning("panel.b")
)

expect_that(tspot.interp(nil = 0, 
                         panel.a = 0, 
                         panel.b = 0, 
                         mito = 0 - tol), 
            gives_warning("mito")
)

})


################################################################################
# tspot.interp warns if given values > 20
test_that("tspot.interp warns if given values > 20", {

expect_that(tspot.interp(nil = 20 + tol, 
                         panel.a = 20, 
                         panel.b = 20, 
                         mito = 20), 
            gives_warning()
)

expect_that(tspot.interp(nil = 20, 
                         panel.a = 20 + tol, 
                         panel.b = 20, 
                         mito = 20), 
            gives_warning()
)

expect_that(tspot.interp(nil = 20, 
                         panel.a = 20, 
                         panel.b = 20 + tol, 
                         mito = 20), 
            gives_warning()
)

expect_that(tspot.interp(nil = 20, 
                         panel.a = 20, 
                         panel.b = 20, 
                         mito = 20 + tol), 
            gives_warning()
)

})



################################################################################
test_that("tspot.cens warns and censors if given values > 20", {

# tspot.cens warns if given values > 20
expect_that(tspot.cens(20 + tol), gives_warning())

# tspot.cens censors if given values > 20
expect_that(tspot.cens(20 + tol), equals(20))

})


################################################################################
# tspot.interp returns a vector of acceptable results
test_that("tspot.interp returns a vector of acceptable results", {

expect_that(tspot.interp(nil = 1:10, 
                         panel.a = 1:10, 
                         panel.b = 1:10, 
                         mito = 1:10), is_a("character"))

})



################################################################################
# tspot.interp tolerates NA values
test_that("tspot.interp tolerates NA values", {

expect_that(
    is.na(tspot.interp(nil = as.numeric(NA), 
                       panel.a = 5, 
                       panel.b = 5, 
                       mito = 20)), is_true()
)

expect_that(
    is.na(tspot.interp(nil = 1, 
                       panel.a = as.numeric(NA), 
                       panel.b = 5, 
                       mito = 20)), is_true()
)

expect_that(
    is.na(tspot.interp(nil = 1, 
                       panel.a = 5, 
                       panel.b = as.numeric(NA), 
                       mito = 20)), is_true()
)

expect_that(
    is.na(tspot.interp(nil = 1, 
                       panel.a = 5, 
                       panel.b = 5, 
                       mito = as.numeric(NA))), is_true()
)

})





