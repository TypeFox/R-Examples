context("indices and positions can be converted to each other")

ind1 <- 3
ind2 <- c(1, 4, 5, 7, 9)
ind3 <- 1:16
ind4 <- c(1, 3, 5, 6, 8, 10, 11, 13, 15)

r1 <- 2
r2 <- 3
r3 <- 4
r4 <- 5

ip1 <- indexToPosition(ind1, 2)
ip2 <- indexToPosition(ind2, 3)
ip3 <- indexToPosition(ind3, 4)
ip4 <- indexToPosition(ind4, 5)

test_that('indexToPosition returns a matrix with 2 columns and length(x) rows', {
    # test that class is a matrix
    expect_is(ip1, 'matrix')
    expect_is(ip2, 'matrix')
    expect_is(ip3, 'matrix')
    expect_is(ip4, 'matrix')
    
    # test there are the proper number of rows
    expect_equal(nrow(ip1), length(ind1))
    expect_equal(nrow(ip2), length(ind2))
    expect_equal(nrow(ip3), length(ind3))
    expect_equal(nrow(ip4), length(ind4))
    
    # test there are two columns
    expect_equal(ncol(ip1), 2)
    expect_equal(ncol(ip2), 2)
    expect_equal(ncol(ip3), 2)
    expect_equal(ncol(ip4), 2)
})

test_that('indexToPosition returns the proper answers', {
    expect_equal(ip1, cbind(row=1, col=2))
    expect_equal(ip2, cbind(row=c(1, 1, 2, 1, 3), col=c(1, 2, 2, 3, 3)))
    expect_equal(ip3, cbind(row=rep(1:4, 4), col=rep(1:4, each=4)))
    expect_equal(ip4, cbind(row=rep(c(1, 3, 5), 3), col=rep(1:3, each=3)))
})


r1 <- 1
r2 <- c(1, 1, 2, 1, 3)
r3 <- rep(1:4, 4)
r4 <- rep(c(1, 3, 5), 3)
c1 <- 2
c2 <- c(1, 2, 2, 3, 3)
c3 <- rep(1:4, each=4)
c4 <- rep(1:3, each=3)

pi1 <- positionToIndex(r1, c1, 2)
pi2 <- positionToIndex(row=r2, col=c2, nrow=3)
pi3 <- positionToIndex(r3, c3, nrow=4)
pi4 <- positionToIndex(r4, c4, nrow=5)

test_that('positionToIndex returns a numeric vector of length row==col', {
    # test class is numeric
    expect_is(pi1, 'numeric')
    expect_is(pi2, 'numeric')
    expect_is(pi3, 'numeric')
    expect_is(pi4, 'numeric')
    
    # test the length equals the number of rows
    expect_equal(length(pi1), length(r1))
    expect_equal(length(pi2), length(r2))
    expect_equal(length(pi3), length(r3))
    expect_equal(length(pi4), length(r4))
    
    # test the length equals the number of columns
    expect_equal(length(pi1), length(c1))
    expect_equal(length(pi2), length(c2))
    expect_equal(length(pi3), length(c3))
    expect_equal(length(pi4), length(c4))
})

test_that('positionToIndex returns the proper answers', {
    expect_equal(pi1, ind1)
    expect_equal(pi2, ind2)
    expect_equal(pi3, ind3)
    expect_equal(pi4, ind4)
})
