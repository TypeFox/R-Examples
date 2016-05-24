
context("DeduImpute")

test_that('deduImpute works for editarrays',{
    E <- editmatrix(c(
        "x1 + x2      == x3",
        "x2           == x4",
        "x5 + x6 + x7 == x8",
        "x3 + x8      == x9",
        "x9 - x10     == x11",
        "x6 >= 0",
        "x7 >= 0"
    ))
    dat <- data.frame(
       x1=c(145,145),
       x2=c(NA,NA),
       x3=c(155,155),
       x4=c(NA,NA),
       x5=c(NA, 86),
       x6=c(NA,NA),
       x7=c(NA,NA),
       x8=c(86,86),
       x9=c(NA,NA),
       x10=c(217,217),
       x11=c(NA,NA)
    )
    v <- deduImpute(E,dat)$corrected
    expect_equal(v$x1,c(145,145))
    expect_equal(v$x2,c(10,10))
    expect_equal(v$x5,c(NA,86))
    expect_equal(v$x6,c(NA,0))
})

test_that('deduImpute handles variables in records not in edits',{
    E <- editmatrix(" x + y == z")
    dat <- data.frame(x=1,y=NA,z=2,v=0)
    v <- deduImpute(E,dat)$corrected
    expect_equal(as.numeric(v[1,]),c(1,1,2,0))
})


context('Deductive imputation with solSpace and imputess')
test_that('solution space works for a simple equality',{
    expect_equal(solSpace(editmatrix("x + y == z"),x=c(x=1,y=NA,z=3))$x0[1],2)
    expect_equal(solSpace(editmatrix("x + y == z"),x=c(x=1,y=NA,z=3))$C[1],0)

})

test_that('solution space works with extra variables in record',{
    expect_equal(solSpace(editmatrix("x + y == z"),x=c(x=1,y=NA,z=3,w=9))$x0[1],2)
    expect_equal(solSpace(editmatrix("x + y == z"),x=c(x=1,y=NA,z=3,u=1,v=NA))$x0[1],2)
})


context('Deductive imputation with deductiveZeros')
test_that('deductiveZeros works for a simple equality',{
    expect_equal(deductiveZeros(editmatrix(c("x + y == z","y>=0")),x=c(x=1,y=NA,z=1)),c(x=FALSE,y=TRUE,z=FALSE))
})

test_that('deductiveZeros works with variables in record not in editmatrix',{
    expect_equal(deductiveZeros(editmatrix(c("x + y == z","y>=0")),x=c(x=1,y=NA,z=1,u=1,v=2)),c(x=FALSE,y=TRUE,z=FALSE,u=FALSE,v=FALSE))
    expect_equal(deductiveZeros(editmatrix(c("x + y == z","y>=0")),x=c(x=1,y=NA,z=1,u=1,v=NA)),c(x=FALSE,y=TRUE,z=FALSE,u=FALSE,v=FALSE))
})

context('The deduImpute method for editset')
test_that('deduImpute.editset works for pure numeric',{
    
    E <- editset(expression(x + y == z))
    x <- data.frame(
        x = NA,
        y = 1,
        z = 1)
    expect_equal(deduImpute(E,x)$corrected$x,0)

})

test_that('deduImpute.editset works for pure categorical',{
    
    E <- editset(expression(
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a' ) B == 'b')
    )
    x <- data.frame(
        A = 'a',
        B = NA)
    expect_equal(deduImpute(E, x)$corrected$B,'b')

})


test_that('deduImpute.editset works for unconnected categorical and numerical',{
    
    E <- editset(expression(
        x + y == z,
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a' ) B == 'b')
    )
    x <- data.frame(
        x = NA,
        y = 1,
        z = 1,
        A = 'a',
        B = NA)
    v <- deduImpute(E,x)
    expect_equal(v$corrected$B,'b')
    expect_equal(v$corrected$x, 0)
    expect_true(v$status$status == 'corrected')
})


test_that('deduImpute.editset works for connected numerical and categorical',{
    E <- editset(expression(
        x + y == z,
        x >= 0,
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a' ) B == 'b',
        if ( B == 'b' ) x > 0
    ))

    x <- data.frame(
        x = NA,
        y = 1,
        z = 1,
        A = 'a',
        B = NA
    )
# this will impute x=0 and B='b',which violates the 
# last edit. Hence, imputation should be reverted.
    v <- deduImpute(E,x) 
    expect_equal(nrow(v$corrections),0)
    expect_equal(v$corrected,x)
})







