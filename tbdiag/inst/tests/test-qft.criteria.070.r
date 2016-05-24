

# Set a tolerance for numeric comparisons
tol <- .Machine$double.eps ^ 0.5


################################################################################
# Standard negatives
test_that("Typical results for negatives are returned as negatives", {

# Low
expect_that(qft.criteria.070(data.frame(nil = 0.01, 
                                        tb = 0.01, 
                                        mito = 10)), 
            matches("Negative")
)

# Borderline
expect_that(qft.criteria.070(data.frame(nil = 0.01, 
                                        tb = 0.70, 
                                        mito = 10)), 
            matches("Negative")
)

# Slightly higher
expect_that(qft.criteria.070(data.frame(nil = 0.50, 
                                        tb = 1.19, 
                                        mito = 10)), 
            matches("Negative")
)

# High, but tb - nil < .25*nil
expect_that(qft.criteria.070(data.frame(nil = 4.00, 
                                        tb = 4.99, 
                                        mito = 10)), 
            matches("Negative")
)

# Nil can be higher than TB
expect_that(qft.criteria.070(data.frame(nil = 0.70, 
                                        tb = 0.01, 
                                        mito = 10)), 
            matches("Negative")
)

# Close to mitogen-indeterminate threshold
expect_that(qft.criteria.070(data.frame(nil = 0.01, 
                                        tb = 0.01, 
                                        mito = 0.51)), 
            matches("Negative")
)

})


################################################################################
# Standard positives
test_that("Typical results for positives are returned as positives", {

# Borderline
expect_that(qft.criteria.070(data.frame(nil = 0.01, 
                                        tb = 0.71, 
                                        mito = 10)), 
            matches("Positive")
)

# Slightly higher
expect_that(qft.criteria.070(data.frame(nil = 0.50, 
                                        tb = 1.20, 
                                        mito = 10)), 
            matches("Positive")
)

# High, tb - nil just over .25*nil
expect_that(qft.criteria.070(data.frame(nil = 4.00, 
                                        tb = 5.00, 
                                        mito = 10)), 
            matches("Positive")
)

# Even if mito - nil < 0.50
expect_that(qft.criteria.070(data.frame(nil = 4.00, 
                                        tb = 5.00, 
                                        mito = 4.25)), 
            matches("Positive")
)

})


################################################################################
# Indeterminate - high nil
test_that("Typical results for high-nil indeterminates are returned as high-nil indeterminates", {

# Nil over 8, everything else standard
expect_that(qft.criteria.070(data.frame(nil = 8.01, 
                                        tb = 0.01, 
                                        mito = 10)), 
            matches("Indeterminate - high nil")
)

# Would be positive if nil weren't over 8
expect_that(qft.criteria.070(data.frame(nil = 8.01, 
                                        tb = 11.00, 
                                        mito = 10)), 
            matches("Indeterminate - high nil")
)

})


################################################################################
# Indeterminate - nil too close to mito
test_that("Typical results for mito~nil indeterminates are returned as mito~nil indeterminates", {

# Lower range
expect_that(qft.criteria.070(data.frame(nil = 0.01, 
                                        tb = 0.01, 
                                        mito = 0.50)), 
            matches("Indeterminate - mitogen too close to nil")
)

# Higher range
expect_that(qft.criteria.070(data.frame(nil = 0.50, 
                                        tb = 0.70, 
                                        mito = 0.99)), 
            matches("Indeterminate - mitogen too close to nil")
)

# Higher still
expect_that(qft.criteria.070(data.frame(nil = 5.00, 
                                        tb = 4.00, 
                                        mito = 5.25)), 
            matches("Indeterminate - mitogen too close to nil")
)

})


################################################################################
# Numeric comparison tolerance checks
test_that("Results are correct in the presence of floating point comparison uncertainty", {

# Nil ~ 8.00
expect_that(qft.criteria.070(data.frame(nil = 8.01, 
                                        tb = 0.01, 
                                        mito = 10)), 
            matches("Indeterminate - high nil")
)

expect_that(qft.criteria.070(data.frame(nil = 8.01 - tol, 
                                        tb = 0.01, 
                                        mito = 10)), 
            matches("Indeterminate - high nil")
)

expect_that(qft.criteria.070(data.frame(nil = 8.01 + tol, 
                                        tb = 0.01, 
                                        mito = 10)), 
            matches("Indeterminate - high nil")
)



# tb - nil ~ 0.70
expect_that(qft.criteria.070(data.frame(nil = 0.13, 
                                        tb = 0.83, 
                                        mito = 10)), 
            matches("Positive")
)

expect_that(qft.criteria.070(data.frame(nil = 0.13 - tol, 
                                        tb = 0.83, 
                                        mito = 10)), 
            matches("Positive")
)

expect_that(qft.criteria.070(data.frame(nil = 0.13 + tol, 
                                        tb = 0.83, 
                                        mito = 10)), 
            matches("Positive")
)



# tb - nil ~ 0.25*nil
expect_that(qft.criteria.070(data.frame(nil = 3.23, 
                                        tb = 3.23 + (3.23 * 0.25), 
                                        mito = 10)), 
            matches("Positive")
)

expect_that(qft.criteria.070(data.frame(nil = 3.23 - tol, 
                                        tb = 3.23 + (3.23 * 0.25), 
                                        mito = 10)), 
            matches("Positive")
)

expect_that(qft.criteria.070(data.frame(nil = 3.23 + tol, 
                                        tb = 3.23 + (3.23 * 0.25), 
                                        mito = 10)), 
            matches("Negative")
)



# mito - nil ~ 0.5
expect_that(qft.criteria.070(data.frame(nil = 0.50, 
                                        tb = 0.01, 
                                        mito = 1.00)), 
            matches("Negative")
)

expect_that(qft.criteria.070(data.frame(nil = 0.50 + tol, 
                                        tb = 0.01, 
                                        mito = 1.00)), 
            matches("Negative")
)

expect_that(qft.criteria.070(data.frame(nil = 0.50 - tol, 
                                        tb = 0.01, 
                                        mito = 1.00)), 
            matches("Negative")
)


})


