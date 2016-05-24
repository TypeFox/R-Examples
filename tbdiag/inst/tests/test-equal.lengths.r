
test_that("equal.lengths accepts equal vecs, errs on unequal, and warns on single vec", {

# equal.lengths is silent when vectors are equal
expect_that(is.null(equal.lengths(1:5, 5:9, 10:14)), is_true())

# Error when they are not equal
expect_that(equal.lengths(1:6, 5:9, 10:14), throws_error())


# Warns when only one vector is passed
expect_that(equal.lengths(1:6), gives_warning())

})
