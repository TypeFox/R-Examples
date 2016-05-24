

context("all_finite")

test_that("all_finite works correctly",{
	expect_true(all_finite(c(1:3)))
	expect_false(all_finite(c(1,2,NA)))
   expect_false(all_finite(c(NA,1,2)))
   expect_false(all_finite(c(NaN,1,2)))
   expect_false(all_finite(c(Inf,1,2)))
   expect_false(all_finite(c(-Inf,1,2)))
	expect_false(all_finite(c(1,2,NA)))
	expect_false(all_finite(c(1,2,NaN)))
	expect_false(all_finite(c(1,2,Inf)))
	expect_false(all_finite(c(1,2,-Inf)))

})




