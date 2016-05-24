test_that("kmers_of", expect_equal(
	kmers_of("hello", 3),
	c("hel", "ell", "llo")))

test_that("vkmers_of", expect_equal(
	vkmers_of(c("hello", "world"), 3),
	list(hello = c("hel", "ell", "llo"), world = c("wor", "orl", "rld"))))
