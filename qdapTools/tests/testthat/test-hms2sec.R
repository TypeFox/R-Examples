context("Checking hms2sec and sec2hms")

test_that("hms2sec and sec2hms produce proper vectors and types",{
	
    a <- c("02:00:03", "04:03:01")
    x <- hms2sec(a)
    b <- c(256, 3456, 56565)
    y <- sec2hms(b)
    
    expect_true(length(x) == length(a))
    expect_true(length(y) == length(b))
    expect_equivalent(x, c(7203, 14581))
    expect_equivalent(y, structure(c(0.00296296296296296, 0.04, 0.6546875), format = "h:m:s", class = "times"))
    
    expect_error(sec2hms(c(NA, a)))
    expect_error(sec2hms(c(NA, b)))
	expect_true(inherits(y, "times"))
})

test_that("hms2sec and sec2hms can be converted between",{
	
    a <- c("02:00:03", "04:03:01")
    x <- hms2sec(a)
    b <- c(256, 3456, 56565)
    y <- sec2hms(b)
    
    expect_true(all(sec2hms(x) == a))
    expect_true(all(hms2sec(y) == b))

})
