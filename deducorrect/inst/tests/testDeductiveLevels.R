
context("Deductive imputation of categorical data")

test_that("deductiveLevels works on simple example",{
    E <- editarray(c(
        "A %in% c('a','b')",
        "B %in% c(TRUE, FALSE)",
        "if ( A == 'a' ) !B"))
    x <- c(A=NA,B=TRUE)
    expect_equal(deductiveLevels(E,x),c(A='b')) 
})

test_that("deductiveLevels works with variables in records/not in editarray",{
    E <- editarray(c(
        "A %in% c('a','b')",
        "B %in% c(TRUE, FALSE)",
        "if ( A == 'a' ) !B"))
    x <- 
    expect_equal(deductiveLevels(E,c(A=NA,B=TRUE,C='c')),c(A='b')) 
    expect_equal(deductiveLevels(E,c(A=NA,B=TRUE,C=NA)),c(A='b')) 

})

test_that("deductiveLevels works on a more complicated example",{

    E <- editarray(c(
        "x1 %in% letters[1:4]",
        "x2 %in% letters[1:3]",
        "x3 %in% letters[1:3]",
        "x4 %in% letters[1:2]",
        "if (x2 == 'c'  & x3 != 'c' & x4 == 'a' ) FALSE",
        "if (x2 != 'a'  & x4 == 'b') FALSE",
        "if (x1 != 'c'  & x2 != 'b' & x3 != 'a') FALSE",
        "if (x1 == 'c'  & x3 != 'a' & x4 == 'a' ) FALSE"
    ))

    x <- c(x1='c',x2='b',x3=NA,x4=NA)
    expect_equal(deductiveLevels(E,x),c(x3='a',x4='a'))

    # another example, partial imputation
    y <- c(x1=NA,x2=NA,x3=NA,x4='b')
    expect_equal(deductiveLevels(E,y),c(x2='a'))
})

test_that("deductiveLevels works on a more complicated example with extra variables",{

    E <- editarray(c(
        "x1 %in% letters[1:4]",
        "x2 %in% letters[1:3]",
        "x3 %in% letters[1:3]",
        "x4 %in% letters[1:2]",
        "if (x2 == 'c'  & x3 != 'c' & x4 == 'a' ) FALSE",
        "if (x2 != 'a'  & x4 == 'b') FALSE",
        "if (x1 != 'c'  & x2 != 'b' & x3 != 'a') FALSE",
        "if (x1 == 'c'  & x3 != 'a' & x4 == 'a' ) FALSE"
    ))

    x <- c(x1='c',x2='b',x3=NA,x4=NA,foo='x')
    expect_equal(deductiveLevels(E,x),c(x3='a',x4='a'))

    # another example, partial imputation
    y <- c(x1=NA,x2=NA,x3=NA,x4='b',bar=NA);
    expect_equal(deductiveLevels(E,y),c(x2='a'))
})



