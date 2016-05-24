context("subMultiple and subVector iterate over vectors well")

theText <- c('Hi Bob & Cooper how is life today', 
             'Anything happening now?', 
             'Sally & Dave are playing with Jess & Julio | with their kids')

subText <- c('Hi Bob and Cooper how is life today', 
             'Anything happening now?', 
             'Sally and Dave are playing with Jess and Julio or with their kids')

multiSub <- subMultiple(theText, pattern=c('&', '\\|'), replacement=c('and', 'or'))
vectSub <- subVector(theText, toSub=c("and"='&', 'or'='\\|'))
vectSubEmpty <- subVector(theText)

#################
## subMultiple
#################

test_that("subMultiple returns the right class", {
    expect_is(multiSub, 'character')
})

test_that("subMultiple returns the right length", {
    expect_equal(length(multiSub), length(theText))
})

test_that("subMultiple gets the right results", {
    expect_identical(multiSub, subText)
})

#################
## subVector
#################

test_that("subVector returns the right class", {
    expect_is(vectSub, 'character')
    expect_is(vectSubEmpty, 'character')
})

test_that("subVector returns the right length", {
    expect_equal(length(vectSub), length(theText))
    expect_equal(length(vectSubEmpty), length(theText))
})

test_that("subVector gets the right results", {
    expect_identical(vectSub, subText)
    expect_identical(vectSubEmpty, theText)
})
