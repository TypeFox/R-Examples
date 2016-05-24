
# Some tests about special preference functionality not tested in test-hasse or test-psel


test_that("Test expression output of preference", {
  
  expect_equal(as.character(as.expression(low(2*a) * -low(b))), 
              "low(2 * a) * -low(b)")
  
})


test_that("Test inherit functions", {
  
  expect_equal(is.preference(empty()), TRUE)
  expect_equal(is.empty_pref(empty()), TRUE)
  expect_equal(is.base_pref(empty()), FALSE)
  
  expect_equal(is.base_pref(low(a) & empty()), TRUE)
  
  expect_equal(is.complex_pref(low(a) & low(b)), TRUE)
})


test_that("Test expression preferences", {
  
  res <- c(30.4, 33.9)
  
  expect_equal(psel(mtcars, low_(expression(hp)) * high(mpg))$mpg, res)
  expect_equal(psel(mtcars, low_(expression(hp * 2)) * high(mpg))$mpg, res)
  expect_equal(psel(mtcars, low(hp) * high_(as.symbol(names(mtcars)[[1]])))$mpg, res)
  
  false <- function(x) -true_(substitute(x))
  expect_equal(as.character(false(cyl == 4)), "-true(cyl == 4)")
  expect_equal(unique(psel(mtcars, false(cyl == 4))$cyl), c(6, 8))
  
})


test_that("Test preferences with df__", {
  
  expect_equal(psel(mtcars, low(df__[[3]]))$disp, 71.1)
  expect_equal(psel(mtcars, high(df__[[4]]))$hp, 335)
  expect_equal(psel(mtcars, true(rownames(df__) == "Fiat 128"))$mpg, 32.4)
  
})


test_that("Test length calculation", {
  
  expect_equal(length(empty()), 0)
  expect_equal(length(layered(a, 1, 2)), 2)
  expect_equal(length(low(a) * -high(b) * empty()), 2)  
  expect_equal(length(true(a) & true(b)), 2)
})


test_that("Test induced layered pref", {
  
  r <- data.frame(A = c(1, 2, 2), B = c(2, 1, 2), C = c(1, 2, 3))
  
  m <- function(a, r) {
    low(psel(r, a, top = nrow(r), show_level = TRUE)[['.level']])
  }

  m_ <- function(a) {
    low(psel(df__, a, top = nrow(df__), show_level = TRUE)$.level)
  }
  
  a <- m(low(A) * low(B), r) & low(C)
  b <- m_(low(A) * low(B)) & low(C)
  
  expect_equal(psel(r, a)$C, 1)
  expect_equal(psel(r, b)$C, 1)
  
  expect_equal(pref.str(a, r), "low(c(1, 1, 2)) & low(C)")
  expect_equal(pref.str(b, r), "low(psel(df__, low(A) * low(B), top = nrow(df__), show_level = TRUE)$.level) & low(C)")
})


test_that("Test evaluations", {
  
  p <- low(f(list(a, c(1,2), c(1,b), list(1,c(2,c)))))
  a <- 1
  b <- 1
  df1 <- data.frame(b=NA, c=NA)
  df2 <- data.frame(c=NA)
  
  expect_equal(as.character(eval.pref(eval.pref(p, df1), df1)), 
              "low(f(list(1, c(1, 2), c(1, b), list(1, c(2, c)))))")
  
  expect_equal(as.character(eval.pref(eval.pref(p, df2), df2)), 
              "low(f(list(1, c(1, 2), c(1, 1), list(1, c(2, c)))))")
  
  
  g <- function(a, ...) low(f(b, ...) + a + sum(...))
  
  expect_equal(as.character(eval.pref(g(1, 2, 3), df1)), 
              "low(f(b, 2, 3) + 1 + 5)")
  
  
                
  f <- function(..., x) prod(...) * x
  g <- function(a, ...) low(f(..., x = b) + (a + sum(...)))
  p <- eval.pref(g(1, 2, 3), df1)
  
  expect_equal(as.character(p), 
               "low(f(2, 3, x = b) + 6)")
  
  expect_equal(psel.indices(data.frame(b=c(1,2)), p), 1)
  
  p <- eval.pref(eval.pref(p, df2), df2)
  
  expect_equal(as.character(p), "low(12)")
  
  expect_equal(psel.indices(data.frame(b=c(1,2)), p), c(1,2))
  
  expect_equal(as.character(true(b == f(2,3,x=4), df1)), 
               "true(b == 24)")
  
  expect_equal(as.character(low(f(b, x = 2 * 2), df1)), 
               "low(f(b, x = 4))")
})
  