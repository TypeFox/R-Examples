library(profilr)
context("frequent_values()")

sample_size <- 10000
set.seed(1)
v <- sample(c(NA, 1:99), size = sample_size, replace = TRUE)
set.seed(1)
w <- sample(c(NA, seq(from = 0, to = 1, length.out = 99)), size = sample_size, replace = TRUE)
set.seed(1)
x <- sample(c(NA, letters), size = sample_size, replace = TRUE)
set.seed(1)
y <- sample(c(NA, TRUE, FALSE), size = sample_size, replace = TRUE)
set.seed(1)
z <- sample(seq(as.Date("1990/1/1"), as.Date("1999/1/1"), "years"), size = sample_size, replace = TRUE)

test_that("most frequent_values work properly ", {
  expect_equal(frequent_values(v, n = 5L, type = "most"), "(1) 25 [125]\n (2) 5 [123]\n (3) 86 [120]\n (4) 40 [119]\n (5) 33 [118]\n")
  expect_equal(frequent_values(w, n = 5L, type = "most"), "(1) 0.244897959183673 [125]\n (2) 0.0408163265306122 [123]\n (3) 0.86734693877551 [120]\n (4) 0.397959183673469 [119]\n (5) 0.326530612244898 [118]\n")
  expect_equal(frequent_values(x, n = 5L, type = "most"), "(1) w [406]\n (2) j [404]\n (3) g [397]\n (4) i [396]\n (5) y [390]\n")
  expect_equal(frequent_values(y, n = 5L, type = "most"), "(1) FALSE [3367]\n (2) TRUE [3310]\n")
  expect_equal(frequent_values(z, n = 5L, type = "most"), "(1) 1999-01-01 [1044]\n (2) 1990-01-01 [1035]\n (3) 1998-01-01 [1024]\n (4) 1994-01-01 [1016]\n (5) 1993-01-01 [1011]\n")
})

test_that("least frequent_values work properly ", {
  expect_equal(frequent_values(v, n = 5L, type = "least"), "(99) 66 [72]\n (97) 46 [76]\n (97) 31 [76]\n (95) 73 [80]\n (95) 60 [80]\n")
  expect_equal(frequent_values(w, n = 5L, type = "least"), "(99) 0.663265306122449 [72]\n (97) 0.459183673469388 [76]\n (97) 0.306122448979592 [76]\n (95) 0.73469387755102 [80]\n (95) 0.602040816326531 [80]\n")
  expect_equal(frequent_values(x, n = 5L, type = "least"), "(26) p [319]\n (25) h [338]\n (24) e [342]\n (22) s [343]\n (22) q [343]\n")
  expect_equal(frequent_values(y, n = 5L, type = "least"), "(2) TRUE [3310]\n (1) FALSE [3367]\n")
  expect_equal(frequent_values(z, n = 5L, type = "least"), "(10) 1996-01-01 [930]\n (9) 1997-01-01 [968]\n (8) 1991-01-01 [977]\n (7) 1995-01-01 [987]\n (6) 1992-01-01 [1008]\n")
})

test_that("frequent_values returns the correct errors", {
  expect_error(frequent_values(x, n = 10L, type = "something"))
  expect_error(frequent_values(x, n = 0L, type = "most"))
  expect_error(frequent_values(x, n = c(1L, 2L), type = "most"))
  expect_error(frequent_values(list(3, 5, 6, 1:10, letters, list(1, 2, 3)), n = 5L, type = "most"))
})


# p <- c(rep(NA, 5), rep(TRUE, 5), rep(FALSE, 5))
# frequent_values(p, n = 5L, type = "most")
