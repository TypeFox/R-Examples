
context("require_a reports errors correctly")

test_that("valid errors are triggered where expected (+ group)", {

	pcall <- "this shows that the expected error was thrown"

	expect_error(require_a(), "missing")

	expect_error(
		require_a(value = 1, pcall = pcall), pcall)
	expect_error(
		require_a(traits = "integer", pcall), pcall)
	expect_error(
		require_a(traits = list(), 1, pcall), pcall)
	expect_error(
		require_a("white-elephant", 1, pcall), pcall)
	expect_error(
		require_a("positive integer", 'string', pcall), pcall)

	expect_error(
		require_a("!null", NULL, pcall), pcall)
	expect_error(
		require_a(c("null", "!unary function"), function (a) a, pcall), pcall)

})

test_that("errors aren't thrown for valid inputs (-)", {

	require_a("positive integer", +1L)
	require_a("whole numeric", 100)
	require_a(c("null", "na", "pairlist"), NULL)
	require_a("matrix", matrix(1:4, 2, 2))

	require_a("integer", 2L, "myfunc(x)")
	require_a("integer", 2L, call('mean', 1,2))

	require_a("!null", 2L, ">_<")
	require_a(c("null", "!unary function"), function (a, b) a+b, "unary()")

})

context("no_match needs more validation")

test_that("no_match displays the call, value, and traits", {

	pcall <- "I am a call"
	value <- "I am the input data"
	traits <- c("positive integer", "matrix")

	expect_error(
		require_a(traits, value, pcall),
		all_patterns(c(
			"positive",
			"integer",
			"matrix",
			value, pcall))
	)
})
