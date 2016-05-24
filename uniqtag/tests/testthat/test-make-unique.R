abc <- c("a", "b", "c")
abcb <- c("a", "b", "c", "b")

test_that("make_unique", expect_equal(
	make_unique(abcb),
	c("a", "b-1", "c", "b-2")))

test_that("make_unique_duplicates", expect_equal(
	make_unique_duplicates(abcb),
	c("a", "b", "c", "b-1")))

test_that("make_unique_all", expect_equal(
	make_unique_all(abcb),
	c("a-1", "b-1", "c-1", "b-2")))

test_that("make_unique_all_or_none 1", expect_equal(
	make_unique_all_or_none(abcb),
	c("a-1", "b-1", "c-1", "b-2")))

test_that("make_unique_all_or_none 2", expect_equal(
	make_unique_all_or_none(abc),
	c("a", "b", "c")))
