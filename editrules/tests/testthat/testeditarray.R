

library(testthat)

context("Editarray")


test_that("2x2 categorical datamodel",{
    dm <- c(
        "g %in% c('m','f')",
        "p %in% c('y','n')")
     expect_equivalent(getArr(editarray(c(dm,"if( p == 'y' )  g != 'm'"))),array(c(F,T,F,T),dim=c(1,4)))
     expect_equivalent(getArr(editarray(c(dm,"if( p %in% c('y') )  g != 'm'"))),array(c(F,T,F,T),dim=c(1,4)))
     expect_equivalent(getArr(editarray(c(dm,"if( p %in% c('y') )  g == 'f'"))),array(c(F,T,F,T),dim=c(1,4)))
})


test_that("2x{TRUE,FALSE} datamodel",{
    dm <- c(
        "g %in% c('m','f')",
        "p %in% c(FALSE,TRUE)")
     expect_equivalent(getArr(editarray(c(dm,"if( p )  g != 'm'"))),array(c(F,T,F,T),dim=c(1,4)))
     expect_equivalent(getArr(editarray(c(dm,"if( g == 'm' ) !p"))),array(c(F,T,F,T),dim=c(1,4)))
     expect_equivalent(getArr(editarray(c(dm,"!p || g=='f'"))),array(c(F,T,F,T),dim=c(1,4)))
})

context("Editarray parsing")
test_that("parse editarray to character and back",{
    edts <- c(
        "g %in% c('m','f')",
        "p %in% c(FALSE,TRUE)",
        "if (p) !g=='m'")
    expect_equivalent(editarray(edts), editarray(as.character(editarray(edts))))
    # cornercase found in version 2.5-1
    edts <- expression(
       A %in% letters[1:3],
       if ( A %in% c('a','b') ) FALSE 
    )
    expect_equivalent(editarray(edts), editarray(as.character(editarray(edts))))
})

test_that("parse editarray to expression and back",{
  edts <- expression(
    g %in% c('m','f'),
    p %in% c(FALSE,TRUE),
    if (p) !g=='m'
    )
  expect_equivalent(editarray(edts), editarray(as.character(editarray(edts))))
})
