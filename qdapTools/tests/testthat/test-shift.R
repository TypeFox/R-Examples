context("Checking shift")

test_that("shift moves left and right",{

    
    x1 <- list(1:10, c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 1L), c(3L, 
        4L, 5L, 6L, 7L, 8L, 9L, 10L, 1L, 2L), c(4L, 5L, 6L, 7L, 8L, 9L, 
        10L, 1L, 2L, 3L), c(5L, 6L, 7L, 8L, 9L, 10L, 1L, 2L, 3L, 4L), 
        c(6L, 7L, 8L, 9L, 10L, 1L, 2L, 3L, 4L, 5L), c(7L, 8L, 9L, 
        10L, 1L, 2L, 3L, 4L, 5L, 6L), c(8L, 9L, 10L, 1L, 2L, 3L, 
        4L, 5L, 6L, 7L), c(9L, 10L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L
        ), c(10L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L))

    x2 <- list(1:10, c(10L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L), c(9L, 
        10L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L), c(8L, 9L, 10L, 1L, 2L, 
        3L, 4L, 5L, 6L, 7L), c(7L, 8L, 9L, 10L, 1L, 2L, 3L, 4L, 5L, 6L
        ), c(6L, 7L, 8L, 9L, 10L, 1L, 2L, 3L, 4L, 5L), c(5L, 6L, 7L, 
        8L, 9L, 10L, 1L, 2L, 3L, 4L), c(4L, 5L, 6L, 7L, 8L, 9L, 10L, 
        1L, 2L, 3L), c(3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 1L, 2L), c(2L, 
        3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 1L))

    expect_equal(lapply(0:9, function(i) shift(1:10, i)), x1)
    expect_equal(lapply(0:9, function(i) shift(1:10, i, "left")), x2)

})

