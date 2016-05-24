
library(dplyr)

# Run all tests for parallel AND for non-parallel mode
for (parallelity in c(FALSE, TRUE)) {
  
  options(rPref.parallel = parallelity)
  
  # Simple tests 
  test_that("Test Preference selection on simple test sets", {
    expect_equal(psel(data.frame(a = c(3,2,1,1,4)), low(a))$a, c(1,1))
    expect_equal(psel.indices(data.frame(a = c(3,3,2,1,1,4)), low(a)), c(4,5))
  })
  
  # Empty preference
  test_that("Test empty preference", {
    expect_equal(psel(data.frame(a = c(3,2,1,1,4)), empty() & high(a))$a, 4)
    expect_equal(psel(data.frame(a = c(3,2,1,1,4)), empty())$a, c(3,2,1,1,4))
    expect_equal(psel(mtcars, low(mpg) & (empty() * low(hp)))$hp, 205)
  })
  
  # Empty dataset
  test_that("Test empty dataset", {
    expect_equal(psel(data.frame(a = 1)[NULL,,drop=FALSE], low(a)), data.frame(a = 1)[NULL,,drop=FALSE])
  })  

  # More tests for psel/psel.indices and the preference constructors
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  test_that("Test Preference selection", {
    expect_equal(psel(mtcars, low(mpg))$mpg, c(10.4, 10.4))
    expect_equal(psel(mtcars, low(mpg) & low(hp))$hp, 205)
    expect_equal(sort(psel(mtcars, high(mpg) * high(hp))$mpg), c(15, 15.8, 17.3, 19.7, 30.4, 32.4, 33.9))
    
    expect_equal(sort(psel(mtcars, true(mpg < 15) * true(am == 0))$mpg), c(10.4, 10.4, 13.3, 14.3, 14.7))
  })
  
  
  test_that("Test Preference selection with indices", {
    expect_equal(sort(psel.indices(mtcars, high(am) * true(vs == 1))), c(3, 18, 19, 20, 26, 28, 32))
  })

  
  # Note that dplyr also uses "between" since v0.3
  test_that("Test base preference macros and prior chains", {
    expect_equal(psel(mtcars, pos(carb, 2))$carb, rep(2, 10))
    expect_equal(psel(mtcars, around(mpg, 20))$mpg, 19.7)
    expect_equal(psel(mtcars, rPref::between(hp, 110, 120))$hp, c(110, 110, 110, 113))
    expect_equal(psel(mtcars, rPref::between(hp, 115, 122))$hp, c(123, 123))
    
    expect_equal(psel(mtcars, -layered(cyl, c(4, 6), 8))$cyl, rep(8, 14))
    expect_equal(rownames(psel(mtcars, true(mpg < 22) & true(cyl == 4) & true(wt < 3 & gear == 4))), "Volvo 142E")
    
  })
  
  test_that("Test if environments are found correctly", {
    test_fun <- function() {
      f <- function(x) -x
      return(low(f(mpg)))
    }
    expect_equal(psel(mtcars, test_fun()), psel(mtcars, -low(mpg)))
  })
  
  
  test_that("Behavior of group_by function from dplyr package", {
    expect_equal(attr(group_by(mtcars[1:5,], cyl), 'indices'), list(2, c(0,1,3), 4))
  })
  
  
  test_that("Grouped preference selection", {
    expect_equal(psel(group_by(mtcars, cyl), low(mpg))$mpg, c(21.4, 17.8, 10.4, 10.4))
    expect_equal(as.data.frame(summarise(psel(group_by(mtcars, cyl), low(mpg) * low(hp)), n())), data.frame(cyl=c(4,6,8),'n()'=c(5,2,2)), check.names=FALSE)
    expect_equal(psel(group_by(mtcars, cc = cyl * carb), true(hp==110) & low(hp))$cc, c(4, 6, 8, 16, 16, 24, 24, 32, 36, 64))
  })
  
  
  # Simple tests of top-K, at_least and toplevel
  test_that("Test TOP-k Preference selection", {
    df <- data.frame(a = c(3,2,1,1,4), b = c(1,1,1,2,2)) # Simple data set
    
    # Check correct indices and level values
    expect_equal(sort(psel.indices(df, low(a), top=5)), 1:5)
    expect_equal(psel(df, low(a), at_least = 2), data.frame(c(1,1), c(1, 2), c(1,1)), check.attributes = FALSE)
    expect_equal(psel.indices(df, low(a), at_least = 3, top = 2), c(3,4))
    expect_equal(psel(df, low(a), top_level = 2)$b, c(1,2,1))
    expect_equal(psel.indices(df, high(a), at_least = 2, top_level = 3, show_level = TRUE)$.indices, c(5,1))
    expect_equal(psel(df, high(a), at_least = 2, top_level = 2, and_connected = FALSE)$.level, c(1,2))
    expect_equal(psel(df, around(a,2), at_least = 2, top = 3, top_level = 2, and_connected = FALSE)$.level, c(1,2,2,2))
    expect_equal(psel(df, around(a,2), at_least = 10)$a, c(2,3,1,1,4))
    expect_equal(psel(df, around(a,2), top_level = 1)$.level, 1)
    expect_equal(psel(df, low(a+b), at_least = 5)$.level, c(1,2,2,3,4))
    
    # Check if show_level works correctly
    expect_equal(psel(df, low(a), show_level = TRUE)$.level, c(1,1))
    expect_equal(ncol(psel(df, low(a))), 2)
    expect_equal(ncol(psel(df, low(a), show_level = TRUE)), 3)
    expect_equal(length(psel.indices(df, low(a))), 2) # ncol is NULL
    expect_equal(length(psel.indices(df, low(a), top_level = 1)), 2)
    expect_equal(ncol(psel.indices(df, low(a), show_level = TRUE)), 2)
    expect_equal(ncol(psel(df, low(a), top = 1)), 3)
    
  })
  
  # Simple tests of grouped top-K, at_least and toplevel
  test_that("Test TOP-k grouped Preference selection", {
    dfg <- group_by(data.frame(a = c(3,2,1,1,4), b = c(1,1,1,2,2)), b) # Simple grouped dataset
    expect_equal(as.data.frame(psel(dfg, low(a), top = 2)), data.frame(c(1,2,1,4), c(1,1,2,2), c(1,2,1,2)), check.attributes = FALSE)
    expect_equal(psel(dfg, low(a), at_least = 2, top = 1)$a, c(1,1))
    expect_equal(psel(dfg, high(a), top_level = 2)$.level, c(1,2,1,2))
    expect_equal(psel.indices(dfg, around(a,2), at_least = 1, top_level = 2), c(2,4))
    expect_equal(psel.indices(dfg, low(b) * high(a), at_least = 2, top_level = 2, show_level = TRUE)$.level, c(1,2,1,2))
  })

  # Top-K Tests on mtcars
  test_that("Test TOP-k Preference selection on mtcars", {
    expect_equal(psel.indices(mtcars, low(mpg + hp), top = 5), order(mtcars$mpg + mtcars$hp)[1:5])
    expect_equal(sort(psel(mtcars, layered(cyl, c(4, 6), 8), top = 4)$cyl), c(4, 6, 6, 6))
    expect_equal(psel(mtcars, low(mpg), top = 5)$mpg, c(10.4, 10.4, 13.3, 14.3, 14.7))
    expect_equal(psel(group_by(mtcars, cyl), low(mpg), top = 3)$mpg, c(21.4, 21.5, 22.8, 17.8, 18.1, 19.2, 10.4, 10.4, 13.3))
    expect_equal(psel(group_by(mtcars, cyl), low(mpg) * high(hp), at_least = 3)$.level, c(1,1,2,1,1,2,2,1,1,1))
  })
  
}  
