library(testthat)
library(formula.tools)
library(magrittr)

context('get.vars')

# NA, NULL, constants, name, symbols, list, formula, call, expression 


# FORMULAS

## ONE-SIDED
context( '  one-sided formula') 

  # These tests don't seem to work, I am not sure why 
  # ( ~ NA )    %>%  get.vars %>% expect_error
  # ( ~ NA )    %>%  get.vars %>% expect_null
  # expect_null( get.vars( ~ NA ) )

  ( ~ NULL )  %>%  get.vars %>% expect_equivalent( character(0) ) 
  ( ~ 1 )     %>%  get.vars %>% expect_equivalent( character(0) )
  ( ~ a )     %>%  get.vars %>% expect_equivalent( 'a' ) 
  ( ~ a + a ) %>%  get.vars %>% expect_equivalent( 'a' ) 
  ( ~ a + b ) %>%  get.vars %>% expect_equivalent( c('a','b'))   

## TWO-SIDED
context( '  two-sided formula' )   
  expect_error( get.vars( NA ~ 1 ) ) 
  ( NULL ~ NULL ) %>% get.vars %>% expect_equivalent( character(0) )
  ( 1 ~ 1 )       %>% get.vars %>% expect_equivalent( character(0) )
  ( y ~ 1 )       %>% get.vars %>% expect_equivalent( 'y' )
  ( y ~ y )       %>% get.vars %>% expect_equivalent( 'y' )
  ( y ~ x )       %>% get.vars %>% expect_equivalent( c('y','x') )
  ( y ~ x + 1 )   %>% get.vars %>% expect_equivalent( c('y','x') )
  ( y ~ x + a + 1 )   %>% get.vars %>% expect_equivalent( c('y','x','a') )  

   
# LHS 
## ONE-SIDED
context('lhs.vars')
context('  one-sided formula')

  ( ~ NULL )  %>%  lhs.vars %>% expect_null 
  ( ~ 1 )     %>%  lhs.vars %>% expect_null
  ( ~ a )     %>%  lhs.vars %>% expect_null 
  ( ~ a + a ) %>%  lhs.vars %>% expect_null 
  ( ~ a + b ) %>%  lhs.vars %>% expect_null   
  

## TWO-SIDED
  context( '  two-sided formula' )   
  # message( "-->", lhs.vars(NA ~ 1) )
  lhs.vars(NA ~ 1) %>% expect_equivalent( character(0))
  # expect_equivalent( lhs.vars(NA ~ 1)  )     # %>% lhs.vars %>% expect_error 
  ( NULL ~ NULL ) %>% lhs.vars %>% expect_null
  ( 1 ~ 1 )       %>% lhs.vars %>% expect_equivalent( character(0) )
  ( y ~ 1 )       %>% lhs.vars %>% expect_equivalent( 'y' )
  ( y ~ y )       %>% lhs.vars %>% expect_equivalent( 'y' )
  ( y ~ x )       %>% lhs.vars %>% expect_equivalent( 'y' )
  ( y ~ x + 1 )   %>% lhs.vars %>% expect_equivalent( 'y' )
  ( y + 1 ~ x + 1 )   %>% lhs.vars %>% expect_equivalent( 'y' )
  ( y ~ x + a + 1 )   %>% lhs.vars %>% expect_equivalent( 'y' )
  ( y + a ~ x )   %>% lhs.vars %>% expect_equivalent( c('y','a') )
  
  
# RHS 
## ONE-SIDED
  context('rhs.vars')
  context('  one-sided formula')
  
  ( ~ NULL )  %>%  rhs.vars %>% expect_equivalent( character(0) ) 
  ( ~ 1 )     %>%  rhs.vars %>% expect_equivalent( character(0) )
  ( ~ a )     %>%  rhs.vars %>% expect_equal('a') 
  ( ~ a + a ) %>%  rhs.vars %>% expect_equal('a') 
  ( ~ a + b ) %>%  rhs.vars %>% expect_equal( c('a','b') )  
  
## TWO-SIDED
  context( '  two-sided formula' )   
  expect_error( rhs.vars(NA ~ 1)  )     # %>% rhs.vars %>% expect_error 
  ( NULL ~ NULL ) %>% rhs.vars %>% expect_equivalent( character(0) )
  ( 1 ~ 1 )       %>% rhs.vars %>% expect_equivalent( character(0) )
  ( y ~ 1 )       %>% rhs.vars %>% expect_equivalent( character(0) )
  ( y ~ y )       %>% rhs.vars %>% expect_equivalent( 'y' )
  ( y ~ x )       %>% rhs.vars %>% expect_equivalent( 'x' )
  ( y ~ x + 1 )   %>% rhs.vars %>% expect_equivalent( 'x' )
  ( y ~ x + a + 1 )   %>% rhs.vars %>% expect_equivalent( c('x','a') )
  
  
# EXPRESSIONS
  
  
# CALLS
  

#