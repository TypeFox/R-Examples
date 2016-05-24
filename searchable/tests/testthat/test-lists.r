library(testthat)
library(searchable)
library(magrittr)

l <- list( ay=1, bee=1:2, cee=1:3, aitch=1:8, aitch=-(1:8) )


context('list-default')
  sl <- searchable(l)

  sl['ee']       %T>% expect_is('list')  %>%  unlist %>% expect_null
  sl['cee']      %>% expect_equivalent( l[3] )# SINGLE HIT
  
  sl[ c('ay','bee','dee') ]  %>% expect_equivalent( l[ c('ay','bee','dee') ] ) 


context('list-regex')
  sl <- searchable(l, 'regex')
  
  sl['ee']  %>% expect_equal( l[2:3] )
  sl['^a'] %>% expect_equal( l[c(1,4,5)])
  
  

context('list-fixed')
  sl <- searchable( l, 'fixed' ) 
  sl['ee']  %>% expect_equal( l[2:3] )


context('list-coll')
  sl <- searchable( l, 'fixed' ) 
  sl['ee'] %>% expect_equal( l[2:3] )
