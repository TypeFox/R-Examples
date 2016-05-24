context("Test interleave")

test_that("interleave works correctly", {
	t1 <- paste0("t1", letters[1:5])
	t2 <- paste0("t2", letters[1:5])
	expect_that(interleave(t1, t2), equals(c("t1a", "t2a", "t1b", "t2b", "t1c", "t2c", "t1d", "t2d", "t1e", "t2e")))
})