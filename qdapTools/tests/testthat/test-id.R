context("Checking id")

test_that("id output is correct length and number of characters",{
    x1 <- id(list(1, 4, 6))
    x2 <- id(matrix(1:10, ncol=1))
    x3 <- id(mtcars)
    x4 <- id(mtcars, TRUE)
    x5 <- id("w")
	x6 <- id(mtcars, prefix="id-")
    
    expect_true(length(x1) == 3 & sum(nchar(x1)) == 3)
    expect_true(length(x2) == 10 & sum(nchar(x2)) == 20)
    expect_true(length(x3) == 32 & sum(nchar(x3)) == 64)
    expect_true(length(x4) == 32 & sum(nchar(x4)) == 128)
    expect_true(length(x5) == 1 & sum(nchar(x5)) == 1)
    expect_true(length(x6) == 32 & sum(nchar(x6)) == 160)
})


