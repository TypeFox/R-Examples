context("Checking pad")

test_that("pad produces vectors of the correct nchar",{
	
    set.seed(10)
    x <- sample(1:10, 10)
    equal.nchar <- function(x){
        a <- nchar(x)
        all(a == a[1])
    }
    
    expect_true(equal.nchar(pad(x)))
    expect_true(equal.nchar(pad(x, sort=FALSE)))
    expect_true(equal.nchar(pad(month.name)))
    expect_true(equal.nchar(pad(month.name, sort=FALSE)))
    expect_true(all(nchar(pad(month.name, padding=20)) == 20))
    expect_false(any(pad(x) == x))
    expect_true(all(as.numeric(pad(x, sort=FALSE)) == x))		
	
})

test_that("pad produces vectors of the correct length",{
	
    set.seed(10)
    x <- sample(1:10, 10)	
	
    expect_true(length(pad(x)) == length(x))
    expect_true(length(pad(x, sort=FALSE)) == length(x))
    expect_true(length(pad(month.name)) == length(month.name))
    expect_true(length(pad(month.name, sort=FALSE)) == length(month.name))
})

test_that("pad produces vectors of the correct mode",{
	
    set.seed(10)
    x <- sample(1:10, 10)
	
    expect_true(mode(pad(x)) == "character")
    expect_true(mode(pad(x, sort=FALSE)) == "character")
    expect_true(mode(pad(month.name)) == "character")
    expect_true(mode(pad(month.name, sort=FALSE)) == "character")
})