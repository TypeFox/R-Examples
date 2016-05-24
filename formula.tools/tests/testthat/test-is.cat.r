library(testthat)
library(formula.tools)
library(magrittr)

context('is.cat')
  # TEST: is.cat
  expect_that( is.cat(letters), is_true() )           # character
  expect_that( is.cat(factor(letters)) , is_true() )  # factor
  expect_that( is.cat(c(TRUE,FALSE)) , is_true() )    # logical
  
  expect_that( is.cat(1:10), is_false() )             # integer
  expect_that( is.cat(1:10/0.5), is_false() )         # numeric
  
  
  # TEST: is.cont
  expect_that( is.cont(letters), is_false() )           # character
  expect_that( is.cont(factor(letters)) , is_false() )  # factor
  expect_that( is.cont(c(TRUE,FALSE)) , is_false() )    # logical
  
  expect_that( is.cont(1:10), is_true() )             # integer
  expect_that( is.cont(1:10/0.5), is_true() )         # numeric
