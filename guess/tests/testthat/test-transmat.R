context("Test transmat")

test_that("transmat works correctly", {
	pre_test_var <- c(1,0,0,1,0,1,0) 
	pst_test_var <- c(1,0,1,1,0,1,1)

	res <- transmat(pre_test_var, pst_test_var)
	cor_ans <- c(x00=2, x01=2, x10=0, x11=3)
	
	expect_that(sapply(res, as.numeric), equals(cor_ans))
})