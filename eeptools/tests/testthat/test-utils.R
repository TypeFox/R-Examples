# Test utils
context("Test defac conversion of factors")

test_that("defac works for all types of factors", {
  a <- as.factor(LETTERS)
  b <- ordered(c(1, 3, '09', 7, 5, "B"))
  expect_is(defac(a), "character")
  expect_is(defac(b), "character")
  a2 <- defac(a)
  b2 <- defac(b)
  expect_identical(levels(a), a2)
  expect_true(all(levels(b) %in% b2))
  expect_identical(length(a), length(a2))
  expect_identical(length(b), length(b2))
})

context("Forcing numerics with makenum")

test_that("makenum works for all types of factors", {
  a <- ordered(c(1, 3, '09', 7, 5))
  a2 <- makenum(a)
  b <- factor(c(1, 3, '09', 7, 5))
  b2 <- makenum(b)
  c <- factor(c(1, 3, '09', 7, 5, "B"))
  c2 <- makenum(c)
  expect_is(a2, "numeric")
  expect_is(b2, "numeric")
  expect_is(c2, "numeric")
  expect_identical(length(a), length(a2))
  expect_identical(length(b), length(b2))
  expect_identical(length(c), length(c2))
  expect_identical(a2, b2)
  expect_identical(c2[6], NA_real_)
})

context("Test that cutoff is numerically accurate")

test_that("cutoff gets the desired result", {
  set.seed(1024)
  a <- rnorm(1000, mean = 0, sd = 1)
  b <- rlnorm(1000, meanlog = 2, sdlog = 1)
  expect_equal(cutoff(a, .05), 0)
  expect_equal(cutoff(a, 0.5), 2)
  expect_equal(cutoff(b, .8), 427)
  
  d <- b
  d[400:500] <- NA
  expect_equal(cutoff(d, 0.2), 131)
  expect_equal(cutoff(d, 0.9, na.rm=FALSE), NA)
  expect_equal(cutoff(d, 0.2, na.rm=FALSE), NA)
  expect_equal(cutoff(d, 0.9, na.rm=TRUE), 648)
  expect_equal(cutoff(d, 0.2, na.rm=TRUE), 131)
  expect_error(cutoff(d, 39))
  expect_error(cutoff(d, -39))
  expect_error(cutoff(d, -0.00039))
})

context("Test the threshold function for numeric accuracy")


test_that("thresh gets the accurate result", {
  set.seed(1024)
  a <- rnorm(1000, mean = 0, sd = 1)
  b <- rlnorm(1000, meanlog = 2, sdlog = 1)
  expect_error(thresh(a, 0))
  expect_equal(thresh(a, 2), 0.5, tol = 0.03, scale = 1)
  expect_equal(thresh(b, 427), 0.8, tol = 0.01)
  
  d <- b
  d[400:500] <- NA
  expect_equal(thresh(d, 131), 0.48, tol = 0.01)
  expect_equal(thresh(d, 648, na.rm=FALSE), NA)
  expect_equal(thresh(d, 131, na.rm=FALSE), NA)
  expect_equal(thresh(d, 600, na.rm=TRUE), 0.92, tol = 0.005)
  expect_equal(thresh(d, 131, na.rm=TRUE), 0.48, tol = 0.01)
  expect_error(thresh(d, 0.39))
  expect_error(thresh(d, -0.39))
  expect_error(thresh(d, -39))
})

context("Test that max_mis works correctly")

test_that("max_mis handles missing data correctly", {
  expect_identical(max(c(7,NA,3,2,0),na.rm=TRUE), max_mis(c(7,NA,3,2,0)))
  max(c(NA,NA,NA,NA),na.rm=TRUE)
  expect_identical(max_mis(c(NA,NA,NA,NA)), NA_real_)
  expect_identical(max_mis(c(NA_real_, NA_real_)), NA_real_)
  expect_identical(max_mis(c()), NA_real_)
  expect_error(max_mis(c("A", "B", "C")))
  expect_error(max_mis(factor("A", "B", "C")))
  expect_error(max_mis(ordered("A", "B", "C")))
})

context("Remove character")

test_that("Remove character works for multiple character type", {
  a <- c(1, 5, 3, 6, "*", 2, 5, "*", "*")
  b <- remove_char(a, "*")
  expect_is(b, "character")
  expect_identical(length(a), length(b))
  expect_equal(length(b[is.na(b)]), 3)
  a <- c(1, 3, 5, "B", "D", ".", ".", ".")
  b <- remove_char(a, ".")
  expect_is(b, "character")
  expect_identical(length(a), length(b))
  expect_equal(length(b[is.na(b)]), 3)
  a <- c(1, 3, 5, "B", "D", "Unk.", "Unk.", "Unk.")
  b <- remove_char(a, "Unk.")
  expect_is(b, "character")
  expect_identical(length(a), length(b))
  expect_equal(length(b[is.na(b)]), 3)
  a <- c(1, 3, 5, "B", "D", "Unk.", "Unk.", "Unk.", NA, NA, NA)
  b <- remove_char(a, "Unk.")
  expect_is(b, "character")
  expect_identical(length(a), length(b))
  expect_equal(length(b[is.na(b)]), 6)
})


context("Leading zero functions as desired")

test_that("Function works for multiple types of inputs", {
  a <- seq(1, 9)
  a2 <- leading_zero(a, digits = 2)
  expect_is(a2, "character")
  expect_true(all(sapply(a2, nchar)==2))
  expect_error(leading_zero(a2, digits = -1))
  expect_error(leading_zero(a2, digits = 0))
  expect_identical(leading_zero(a, digits = -1), leading_zero(a, digits = 0))
  
  a <- seq(9, 25)
  a2 <- leading_zero(a, digits = 3)
  expect_is(a2, "character")
  expect_true(all(sapply(a2, nchar)==3))
  a2 <- leading_zero(a, digits = 1)
  expect_false(all(sapply(a2, nchar)==1))
  
  expect_error(leading_zero(c("A", "B", "C", digits = 2)))
  a <- c(-5000, -50, -5, -.01, 0, 0.1, 4, 40, 400, 4000)
  a2 <- leading_zero(a, digits = 3)
  expect_identical(a2, c("-5000", "-050", "-005", "0000", "0000", "0000", 
                         "0004", "0040", 
                         "0400", "4000"))
})

context("Test decomma")

a <- c("12,124", "21,131", "A,b")
b <- c("12124", "21131", "Ab")

c <- a[1:2]
d <- as.numeric(b[1:2])

test_that("decomma returns the right class", {
  expect_that(decomma(c), equals(d))
  expect_that(decomma(a), gives_warning())
  expect_that(decomma(a), is_a("numeric"))
  expect_that(decomma(b), is_a("numeric"))
  expect_that(decomma(c), is_a("numeric"))
  expect_that(decomma(d), is_a("numeric"))
})

n <- c(NA, NA, NA, "7,102", "27,125", "23,325,22", "Ab")

test_that("decomma handles NAs properly", {
  expect_that(length(decomma(n)[!is.na(decomma(n))]), equals(3))
  expect_that(decomma(n)[6], equals(2332522))
})


context("nth max")

test_that("Numeric accuracy", {
  a <- c(1:20, 20:1)
  b <- c(LETTERS[c(2, 9, 3, 12, 1)])
  z <- c(121253125, 12401892377905, 31221, 12, 45, -2145125, -123, 0)
  f <- c(10, 10, 10, 10, 9, 9, 10.0001, 10.0001)
  expect_error(nth_max(a, 100), "index .* outside bounds")
  expect_error(nth_max(a, -1), "index .* outside bounds")
  expect_equal(nth_max(a), 20)
  expect_equal(nth_max(b, 3), "C")
  expect_equal(nth_max(z), 12401892377905)
  expect_equal(nth_max(f), 10.0001)
})
