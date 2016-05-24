library(testthat)
library(searchable)
library(magrittr)

v <- c( ay=1, bee=2, cee=3 )


context('vector-default')
  sv <- searchable(v)
  
  sv['cee']    %>% expect_equivalent(3)
  sv['ee']     %>% is.na %>% expect_true
  sv[ c('ay','bee','dee') ]  %>% expect_equivalent( c(1,2,NA) ) 


context('vector-regex')
  sv <- v %>% searchable('regex')
  
  sv@type %>% expect_equal('regex')

  sv['ay'] %>% expect_equal( v['ay'] )
  sv['ee'] %>% expect_equal( v[2:3] )
  sv['e$'] %>% expect_equal( v[2:3] )
  

context('vector-fixed')
  sv <- searchable(v, 'fixed' )

  sv@type %>% expect_equal('fixed')
  
  sv['cee']    %>% expect_equivalent(3)
  sv['ee']     %>% expect_equivalent(2:3)

context('vector-coll')


context('ignore-case')

  sv <- v %>% searchable('regex', case_insensitive = FALSE ) 
  sv@options$case_insensitive %>% expect_false

  sv <- sv %>% ignore.case 
  sv@options$case_insensitive %>% expect_true
  
  sv <- sv %>% use.case 
  sv@options$case_insensitive %>% expect_false
  
  