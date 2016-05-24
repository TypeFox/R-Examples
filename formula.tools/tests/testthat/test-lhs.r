library(testthat)
library(formula.tools)
library(magrittr)


# FORMULA
context( 'formula')

## LHS
context('  lhs')

  lhs( NULL ~ . ) %>% expect_null
  lhs( NA ~ . )   %>% expect_equal(NA) 
  lhs( . ~ . )    %>% expect_equal( as.name('.') )
  lhs( 1 ~ . )    %>% expect_equal(1)
  lhs( a ~ . )    %>% expect_equal( as.name('a') )
  lhs( a+1 ~ .  ) %>% expect_equivalent( expression(a+1)[[1]] )
  lhs( a+b ~ .  ) %>% expect_equivalent( expression(a+b)[[1]] )
  lhs( ~ .  )     %>% expect_null
  
## RHS
context('  rhs')

  rhs( . ~ NULL  ) %>% expect_null
  rhs( . ~ NA )   %>% expect_equal(NA) 
  rhs( . ~ . )    %>% expect_equal( as.name('.') )
  rhs( . ~ 1 )    %>% expect_equal(1)
  rhs( . ~ a  )   %>% expect_equal( as.name('a') )
  rhs( . ~ a+1 )  %>% expect_equivalent( expression(a+1)[[1]] )
  rhs( . ~ a+b )  %>% expect_equivalent( expression(a+b)[[1]] )
  rhs( ~ .  )     %>% expect_equal( as.name('.') )


# FORMULA
context( 'call')
  
  ## LHS
  context('  lhs')
  
  expression(NULL + .)[[1]]  %>% lhs %>% expect_null
  expression( NA + . )[[1]]  %>% lhs %>% expect_equal(NA) 
  expression( . + . )[[1]]   %>% lhs %>% expect_equal( as.name('.') )
  expression( 1 + . )[[1]]   %>% lhs %>% expect_equal(1)
  expression( a + . )[[1]]   %>% lhs %>% expect_equal( as.name('a') )
  expression( a+1 + . )[[1]] %>% lhs %>% expect_equivalent( expression(a+1)[[1]] )
  expression( a+b + . )[[1]] %>% lhs %>% expect_equivalent( expression(a+b)[[1]] )
  expression( ~ . )[[1]]     %>% lhs %>% expect_null
  
  ## RHS
  context('  rhs')
  
  expression( . / NULL )[[1]] %>% rhs %>% expect_null
  expression( . / NA )[[1]]   %>% rhs %>% expect_equal(NA) 
  expression( . / . )[[1]]    %>% rhs %>% expect_equal( as.name('.') )
  expression( . / 1 )[[1]]    %>% rhs %>% expect_equal(1)
  expression( . / a  )[[1]]   %>% rhs %>% expect_equal( as.name('a') )
  expression( ./a + 1 )[[1]]  %>% rhs %>% expect_equivalent( 1 )
  expression( ./a + b )[[1]]  %>% rhs %>% expect_equivalent( as.name('b') )
  expression( ~ . )[[1]]      %>% rhs %>% expect_equal( as.name('.') )
