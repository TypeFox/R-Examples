context("Flag usage")

test_that("Dictionary switcher works",{
	test <- list("if length(out) == 2", "if 2 == 4L", "if 1 == 3L", "if 3", "if finally")
	swc <- makeFunSwitcher(test)
	
	expect_equal(swc(c("meh", "badSort"), numOut = 2), 1)
	expect_equal(swc(c("goodSort", "4")), 2)
	expect_equal(swc(c("3", "badSort")), 3)
	expect_equal(swc(c("meh", "badSort", "do this too")), 4)
	expect_equal(swc(c("meh", "badSort")), 5)
})
