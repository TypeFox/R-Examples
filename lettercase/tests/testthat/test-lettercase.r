library(lettercase)
library(testthat)


# INITIALIZATION
  context( 'LetterCase')
#   fooBar <- c( 'foo', 'BAR' )
#   
#   foo_bar <- LetterCase( fooBAR )
#   
#   
#   expect_is( foo_bar, 'LetterCase' )
#   expect_equal( length(foo_bar), 2 )
#   
#   expect_that( foo_bar@.Data, is_equivalent_to( c('foo', 'BAR' ) ) )
# 
# 
#   
# # NAMED ELEMENTS
# 
#   context( 'Named elements')
#   fooBar <- c( a='foo', b='BAR' )
# 
#   foo_bar <- LetterCase( fooBar )
# 
#   expect_is( foo_bar, 'LetterCase' )
#   expect_equal( length(foo_bar), 2 )
#   
#   expect_that( foo_bar@.Data, is_equivalent_to( c('foo', 'BAR' ) ) )
# 
#   
# # ACCESSORS
# 
# 
# # SUBCLASS
#   context( 'LetterCase subclass')
#   
#   SubLetterCase <- setClass( 
#     'SubLetterCase', contains='LetterCase'
#     , prototype = prototype(.init=I) 
#   ) 
# 
#   foo_bar <- SubLetterCase( fooBAR )
# 
#   expect_is( foo_bar, 'SubLetterCase')
#   expect_equal( length(foo_bar), 2 )
# 
#   expect_that( foo_bar@.Data, is_equivalent_to( c('foo', 'BAR' ) ) )
# 
# 
# # UPPERCASE
# 
#   context( 'UpperCase' )
#   UpperCase <- setClass( 
#     'UpperCase', contains='LetterCase'
#     # , prototype = prototype( .init=toupper )
#   )
# 
#   UpperCase(letters)
