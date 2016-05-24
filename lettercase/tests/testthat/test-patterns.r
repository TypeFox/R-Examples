library(testthat)
library(stringr)

context('patterns')
  whitespace <- c( " ", "  ", "\t", "\t\t")
  for( ws in whitespace ) 
    expect_true( grepl( pattern_whitespace, ws ) )

  whitespace_like <- c( "-", "_" )
  for( ws in whitespace_like ) 
    expect_true( grepl( pattern_whitespace_like, ws ) )