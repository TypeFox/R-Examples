library(testthat) 
library(magrittr)
library(operator.tools)

context("fun2name")

  fun2name( ls ) %>% expect_equal( 'ls' )
  fun2name( `%in%` ) %>% expect_equal( '%in%' )
  
context("name2fun")
  name2fun( 'ls' ) %>% expect_identical( ls )
  name2fun( '%in%' ) %>% expect_identical( `%in%` )
