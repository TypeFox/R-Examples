
context ("ensure that custom trait tests work as expected")

test_that("finite works", {

	f <- trait_tests$finite
	expect_false(f(+Inf))
	expect_false(f(-Inf))
	expect_true(f(10))
})


test_that("infinite works", {

	f <- trait_tests$infinite
	expect_true(f(+Inf))
	expect_true(f(-Inf))
	expect_false(f(10))
})


test_that("na works", {

	f <- trait_tests$na
	expect_true(f(NA))
	expect_false(f(10))
})


test_that("nan works", {

	f <- trait_tests$nan
	expect_true(f(NaN))
	expect_false(f(100))
})


test_that("object works", {

	f <- trait_tests$object
	expect_true( f(structure(list(a = 1), class = "me")) )
	expect_false( f(list(a = 1)) )
})


test_that("false works", {

	f <- trait_tests$false
	expect_true(f(FALSE))
	expect_false(f(TRUE))
	expect_false(f(NA))
})


test_that("positive works", {

	f <- trait_tests$positive
	expect_false(f(0))
	expect_true(f(1))
	expect_false(f(-10))
})


test_that("nonnegative works", {

	f <- trait_tests$nonnegative
	expect_true(f(0))
	expect_true(f(1))
	expect_false(f(-10))
})

test_that("whole works", {

	f <- trait_tests$whole
	expect_true(f(1))
	expect_true(f(0L))
	expect_false(f(1.000001))
})

test_that("named works", {

	f <- trait_tests$named
	expect_true( f(list(a = 1, b = 2)) )
	expect_false( f(list(a = 1, 2)) )
	expect_true( f(c(a = 1, b = 2)) )
	expect_false( f(c(a = 1, 2)) )
	expect_false( f(10) )
})

test_that("functionable works", {

	f <- trait_tests$functionable

	expect_true( f("a") )
	expect_true( f(as.symbol("a")) )
	expect_true( f(mean) )
	expect_false( f(1) )
	expect_false( f(list(1, 2, 3)) )

})

test_that("boolean works", {

	f <- trait_tests$boolean
	expect_false(f(NA))
	expect_true(f(TRUE))
	expect_true(f(FALSE))
})


test_that("string works", {

	f <- trait_tests$string
	expect_true(f("this"))
	expect_false(f(c("a", "b")))
	expect_false(f(+100))
})


test_that("listy works", {

	f <- trait_tests$listy
	expect_true( f(pairlist(1)) )
	expect_true( f(list(1)) )
	expect_true(f(1))
})


test_that("length_zero works", {

	f <- trait_tests$length_zero
	expect_true(f(NULL))
	expect_false(f(1))
	expect_false(f(1:2))
})


test_that("length_one works", {

	f <- trait_tests$length_one
	expect_false(f(NULL))
	expect_true(f(1))
	expect_false(f(1:2))
})


test_that("nullary works", {

	f <- trait_tests$nullary
	expect_true(f(function () x))
	expect_true(f(function (...) x))
	expect_false(f(function (x, y) x))

})

test_that("unary works", {

	f <- trait_tests$unary
	expect_true(f(function (x) x))
	expect_true(f(function (...) x))
	expect_false(f(function (x, y) x))

})

test_that("variadic works", {

	f <- trait_tests$variadic
	expect_true(f( function (...) x))
	expect_false(f(function (x, y) x))
	expect_true(f(function (x, ...) x))

})


