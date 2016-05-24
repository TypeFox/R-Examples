context("Test Standard Guessing Correction")

test_that("Standard Guessing Correction Happens Correctly", {
	pre_test <- data.frame(item1=c(1,0,0,1,0), item2=c(1,NA,0,1,0)) 
	pst_test <-  pre_test + cbind(c(0,1,1,0,0), c(0,1,0,0,1))
	lucky <- rep(.25, 2)
	# Adjusted Effect
	res <- stndcor(pre_test, pst_test, lucky)
	expect_that(res, is_a("list"))
})