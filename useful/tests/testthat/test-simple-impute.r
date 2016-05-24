context('Checking that simple impute does as intended.')

# Create fake data with missing values
theDF <- data.frame(A=1:10, B=1:10, C=1:10)
theDF[c(1, 4, 6), c(1)] <- NA
theDF[c(3, 4, 8), c(3)] <- NA

test_that('Median works', {
    expect_true(is.na(median(theDF$A)))
    expect_equal(median(theDF$A, na.rm=TRUE), expected=7)
    expect_true(!is.na(median(theDF$B)))
    expect_equal(median(theDF$B, na.rm=TRUE), expected=5.5)
    expect_true(is.na(median(theDF$C)))
    expect_equal(median(theDF$C, na.rm=TRUE), expected=6)
})

test_that('simple.impute.default properly imputes the median', {
    expect_equal(simple.impute.default(theDF$A), c(7, 2, 3, 7, 5, 7, 7, 8, 9, 10))
    expect_equal(simple.impute.default(theDF$B), 1:10)
    expect_equal(simple.impute.default(theDF$C), c(1, 2, 6, 6, 5, 6, 7, 6, 9, 10))
    
    expect_equal(simple.impute(theDF$A), c(7, 2, 3, 7, 5, 7, 7, 8, 9, 10))
    expect_equal(simple.impute(theDF$B), 1:10)
    expect_equal(simple.impute(theDF$C), c(1, 2, 6, 6, 5, 6, 7, 6, 9, 10))
    
    expect_equal(length(simple.impute(theDF$A)), 10)
    expect_equal(length(simple.impute(theDF$B)), 10)
    expect_equal(length(simple.impute(theDF$C)), 10)
})

test_that('simple.impute.default properly imputes the mean', {
    expect_output(all.equal(simple.impute.default(theDF$A, mean), c(6.285714, 2.000000, 3.000000, 6.285714, 5.000000, 6.285714, 7.000000, 8.000000, 9.000000, 10.000000)),
                  "Mean relative difference: 4.545455e-08")
    expect_equal(simple.impute.default(theDF$B), 1:10)
    expect_output(all.equal(simple.impute.default(theDF$C, mean), c(1.0, 2.0, 5.714286, 5.714286, 5.0, 6.0, 7.0, 5.714286, 9.0, 10.0)), 
                  "Mean relative difference: 5e-08")
    
    expect_output(all.equal(simple.impute(theDF$A, mean), c(6.285714, 2.000000, 3.000000, 6.285714, 5.000000, 6.285714, 7.000000, 8.000000, 9.000000, 10.000000)),
                  "Mean relative difference: 4.545455e-08")
    expect_equal(simple.impute(theDF$B), 1:10)
    expect_output(all.equal(simple.impute(theDF$C, mean), c(1.0, 2.0, 5.714286, 5.714286, 5.0, 6.0, 7.0, 5.714286, 9.0, 10.0)), 
                  "Mean relative difference: 5e-08")
    
    expect_equal(length(simple.impute(theDF$A, mean)), 10)
    expect_equal(length(simple.impute(theDF$B, mean)), 10)
    expect_equal(length(simple.impute(theDF$C, mean)), 10)
})

test_that('simple.impute.default properly imputes a constant', {
    expect_equal(simple.impute.default(theDF$A, constant(4)), c(4, 2, 3, 4, 5, 4, 7, 8, 9, 10))
    expect_equal(simple.impute.default(theDF$B, constant(4)), 1:10)
    expect_equal(simple.impute.default(theDF$C, constant(4)), c(1, 2, 4, 4, 5, 6, 7, 4, 9, 10))
    
    expect_equal(simple.impute(theDF$A, constant(4)), c(4, 2, 3, 4, 5, 4, 7, 8, 9, 10))
    expect_equal(simple.impute(theDF$B, constant(4)), 1:10)
    expect_equal(simple.impute(theDF$C, constant(4)), c(1, 2, 4, 4, 5, 6, 7, 4, 9, 10))
    
    expect_equal(length(simple.impute(theDF$A, constant(4))), 10)
    expect_equal(length(simple.impute(theDF$B, constant(4))), 10)
    expect_equal(length(simple.impute(theDF$C, constant(4))), 10)
    
    expect_equal(simple.impute(theDF$A, constant()), c(1, 2, 3, 1, 5, 1, 7, 8, 9, 10))
    expect_equal(simple.impute(theDF$B, constant()), 1:10)
    expect_equal(simple.impute(theDF$C, constant()), c(1, 2, 1, 1, 5, 6, 7, 1, 9, 10))
})

test_that('simple.impute.data.frame works as expected', {
    expect_equal(dim(simple.impute(theDF)), c(10, 3))
    expect_equal(simple.impute(theDF), 
                 data.frame(A=simple.impute(theDF$A), 
                            B=simple.impute(theDF$B), 
                            C=simple.impute(theDF$C)
                 )
    )
    expect_equal(simple.impute(theDF, mean), 
                 data.frame(A=simple.impute(theDF$A, mean), 
                            B=simple.impute(theDF$B, mean), 
                            C=simple.impute(theDF$C, mean)
                 )
    )
    expect_equal(simple.impute(theDF, constant(4)), 
                 data.frame(A=simple.impute(theDF$A, constant(4)), 
                            B=simple.impute(theDF$B, constant(4)), 
                            C=simple.impute(theDF$C, constant(4))
                 )
    )
})

test_that('constant returns a function', {
    expect_is(constant(4), 'function')
})

test_that('constant(n)(x) results in n', {
    expect_equal(constant(4)(1:10), 4)
    expect_equal(constant(1)(1:10), 1)
    expect_equal(constant()(1:10), 1)
    expect_equal(constant(5)(1), 5)
})