


################################################################################
# Standard negatives
test_that("Typical results for negatives are returned as negatives", {

# Low
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 1, 
                                                  panel.a = 1, 
                                                  panel.b = 1, 
                                                  mito = 20)), 
            matches("Negative")
)

# Just below borderline - Panel A
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 1, 
                           panel.a = 5, 
                           panel.b = 1, 
                           mito = 20)), 
            matches("Negative")
)

# Just below borderline - Panel B
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 1, 
                           panel.a = 1, 
                           panel.b = 5, 
                           mito = 20)), 
            matches("Negative")
)


# Slightly higher
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                           panel.a = 9, 
                           panel.b = 9, 
                           mito = 20)), 
            matches("Negative")
)


# Nil can be higher than TB
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 8, 
                           panel.a = 1, 
                           panel.b = 1, 
                           mito = 20)), 
            matches("Negative")
)

})


################################################################################
# Standard positives
test_that("Typical results for positives are returned as positives", {

# Low edge of positive - Panel A
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 1, 
                                       panel.a = 9, 
                                       panel.b = 6, 
                                       mito = 20)), 
            matches("Positive")
)

# Low edge of positive - Panel B
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 1, 
                                       panel.a = 6, 
                                       panel.b = 9, 
                                       mito = 20)), 
            matches("Positive")
)

# Slightly higher - Panel A
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 16, 
                                       panel.b = 6, 
                                       mito = 20)), 
            matches("Positive")
)

# Slightly higher - Panel B
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 6, 
                                       panel.b = 16, 
                                       mito = 20)), 
            matches("Positive")
)


# Even when mito < 20
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 6, 
                                       panel.b = 16, 
                                       mito = 10)), 
            matches("Positive")
)




})



################################################################################
# Borderline

# 5 spots - Panel A
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 10, 
                                       panel.b = 1, 
                                       mito = 20)), 
            matches("Borderline")
)


# 6 spots - Panel A
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 11, 
                                       panel.b = 1, 
                                       mito = 20)), 
            matches("Borderline")
)


# 7 spots - Panel A
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 12, 
                                       panel.b = 1, 
                                       mito = 20)), 
            matches("Borderline")
)


# 5 spots - Panel B
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 1, 
                                       panel.b = 10, 
                                       mito = 20)), 
            matches("Borderline")
)


# 6 spots - Panel B
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 1, 
                                       panel.b = 11, 
                                       mito = 20)), 
            matches("Borderline")
)


# 7 spots - Panel B
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 1, 
                                       panel.b = 12, 
                                       mito = 20)), 
            matches("Borderline")
)


# Even when mito < 20
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 12, 
                                       panel.b = 1, 
                                       mito = 10)), 
            matches("Borderline")
)

expect_that(tspot.criteria.oxford.usa(data.frame(nil = 5, 
                                       panel.a = 1, 
                                       panel.b = 12, 
                                       mito = 10)), 
            matches("Borderline")
)





################################################################################
# Invalid - high nil
test_that("Typical results for high-nil invalids are returned as high-nil invalids", {

# Nil over 10, everything else standard
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 11, 
                                       panel.a = 1, 
                                       panel.b = 1, 
                                       mito = 20)), 
            matches("Invalid - high nil")
)

# Takes precedence over low mitogen
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 11, 
                                       panel.a = 1, 
                                       panel.b = 1, 
                                       mito = 10)), 
            matches("Invalid - high nil")
)

})


################################################################################
# Invalid - low mito
test_that("Typical results for low-mito invalids are returned as low-mito invalids", {

# Very low
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 1, 
                                      panel.a = 1, 
                                      panel.b = 1, 
                                      mito = 10)), 
            matches("Invalid - low mitogen")
)

# At the threshold
expect_that(tspot.criteria.oxford.usa(data.frame(nil = 1, 
                                      panel.a = 1, 
                                      panel.b = 1, 
                                      mito = 19)), 
            matches("Invalid - low mitogen")
)

})





################################################################################
# Results match known data
load(file.path("..", "..", "data", "test.tspots.rdata"))

# Compute results
fun.result <- trim.output(tspot.criteria.oxford.usa(test.tspots), "terse")

# Compare to lab results
test_that("tspot.criteria.oxford.usa results exactly match original lab results", {

expect_that(all.equal(test.tspots$lab.result, fun.result), is_true())

})
