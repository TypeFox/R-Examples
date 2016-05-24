context("Checking start_end")

test_that("start_end gets locations of start/end places for the ones in a binary vector",{
	
    set.seed(10); x <- sample(0:1, 50, TRUE, c(.35, .65))
    expected <- structure(list(start = c(1L, 5L, 12L, 22L, 24L, 28L, 30L, 35L, 
        40L, 45L), end = c(3L, 10L, 19L, 22L, 25L, 28L, 33L, 35L, 43L, 
        49L)), .Names = c("start", "end"), row.names = c(NA, -10L), class = "data.frame")
    
    expect_equal(expected, start_end(x))
	
})
