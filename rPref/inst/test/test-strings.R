

test_that("Test string output of preference", {

  expect_output(show.pref(empty()), '[Preference] (empty)', fixed = TRUE)
  
  expect_output(show.pref(low(a) * low(b)), 
                '[Preference] low(a) * low(b)', fixed = TRUE)
  
  expect_output(show.pref(low(a) * (low(b) & high(c))), 
                '[Preference] low(a) * (low(b) & high(c))', fixed = TRUE)
  
  expect_output(show.pref((low(a) * low(b)) & high(c)), 
                '[Preference] (low(a) * low(b)) & high(c)', fixed = TRUE)
  
})


test_that("Test SQL output", {
  
  expect_identical(show.query(empty()), "")
  
  expect_identical(show.query(low(a1 * a2) * high(b) * (-high(c) & true(d)), dialect = "EXASOL"),
                   "PREFERRING LOW (a1 * a2) PLUS HIGH b PLUS (INVERSE (HIGH c) PRIOR TO d)")
  
  expect_identical(show.query(low(a1 * a2) * high(b) * (-high(c) & true(d)), dialect = "Preference SQL"),
                   "PREFERRING (a1 * a2) LOWEST AND b HIGHEST AND ((c HIGHEST) DUAL PRIOR TO d = TRUE)")
  
})


test_that("Test string output of preferences on a given data set", {
  
  f <- function(x) (2*x)
  y <- 1
  
  expect_output(show.pref(low(wt) * low(hp) * low(f(cyl + y)), df = mtcars), 
                '[Preference] low(wt) * low(hp) * low(f(cyl + 1))', fixed = TRUE)
  
  expect_output(show.pref(eval.pref(low(wt) * low(hp) * true(f(cyl + y) > y + y), df = mtcars)), 
                '[Preference] low(wt) * low(hp) * true(f(cyl + 1) > 2)', fixed = TRUE)
  
  expect_identical(pref.str((low(wt) * low(hp)) & reverse(high(y + f(cyl))), df = mtcars), 
                  '(low(wt) * low(hp)) & -high(1 + f(cyl))')
  
  expect_identical(pref.str((low(wt) * low(hp)) & reverse(high(y + f(cyl))), df = mtcars), 
                   '(low(wt) * low(hp)) & -high(1 + f(cyl))')
  
  expect_identical(as.character(eval.pref(-high(f(y) + f(cyl)), df = mtcars)), 
                   '-high(2 + f(cyl))')
  
  expect_identical(show.query((low(wt) * low(hp)) & high(cyl + f(wt + y)), df = mtcars),
                   "PREFERRING (LOW wt PLUS LOW hp) PRIOR TO HIGH (cyl + f(wt + 1))")
  
  expect_identical(show.query((low(wt) * low(hp)) | high(cyl + f(wt + y)), df = mtcars, dialect = "Psql"),
                  "PREFERRING (wt LOWEST AND hp LOWEST) INTERSECT WITH (cyl + f(wt + 1)) HIGHEST")
})
