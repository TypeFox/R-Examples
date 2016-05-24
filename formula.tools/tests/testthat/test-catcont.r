library(testthat)
library(formula.tools)
library(magrittr)


context( "is.cat")

  is.cat(letters)          %>% expect_true
  is.cat(factor(letters))  %>% expect_true
  is.cat(TRUE)             %>% expect_true 
  is.cat(FALSE)            %>% expect_true 
  is.cat(1:10)             %>% expect_false 
  is.cat(rnorm(10))        %>% expect_false 
  # is.cat( now() )          %>% expect_false

context('is.cont')  
  is.cont(letters)         %>% expect_false
  is.cont(factor(letters)) %>% expect_false  
  is.cont(TRUE)            %>% expect_false 
  is.cont(FALSE)           %>% expect_false
  is.cont(1:10)            %>% expect_true
  is.cont(rnorm(10))       %>% expect_true  
  
context('which.cat')
  which.cat(iris)               %>% expect_equal( c( Species= 5) )
  which.cat( iris, names=TRUE ) %>% expect_equal( "Species")
  
  
context('which.cont')
  which.cont( iris )             %>% expect_equivalent( 1:4 )  
  which.cont( iris, names=TRUE ) %>% expect_equal( names(iris)[1:4] )
  
