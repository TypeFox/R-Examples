library(testthat)
library(operator.tools)

context( 'Function: is.operator' )
for( op in operators() ) 
  expect_that( 
    is.operator( as.name(op)), 
    is_true(), 
    paste( op, "is not an operator." ) 
  )


context( 'Function: can.operator' )
for( op in operators() ) 
  expect_that( 
    can.operator(op), 
    is_true(), 
    paste( op, "cannot be coercised into an operator." ) 
  )



context('Function: as.operator')
for( op in operators() )
  expect_that( 
    is.operator( as.operator(op) ),  
    is_true() ,
    paste( op, 'was not successfully coerced to an operator.' )
  )


