

context("Detecting blocks")


test_that("editmatrix separates in blocks",{
    E  <- editmatrix( c( 
        "x1 + x2 == x3"
      , "x3 + x4 == x5"
      , "x5 + x6 == x7"
      , "y1 + y2 == y3"
      , "z1 + z2 == z3")
    )
    expect_equal(length(blocks(E)),3)
    expect_true(all(getAb(blocks(E)[[2]])==c(1,1,-1,0)))    
})

test_that("editarray separates in blocks",{
    E <- editarray(c(
        "x %in% c('a','b','c')",
        "y %in% c('d','e')",
        "z %in% c('f','g')",
        "u %in% c('w','t')",
        "if ( x == 'a') y != 'd'",
        "if ( z == 'f') u != 'w'"))
    expect_equal(length(blocks(E)),2)
    expect_true(all(getArr(blocks(E)[[1]]) == c(TRUE,FALSE,FALSE,TRUE,FALSE)))
})

test_that("editset separates in blocks",{
    E <- editset(expression(
        if ( x > 0 ) y > 0,
        x + y >= z,
        A %in% letters[1:2],
        B %in% letters[2:3],
        if ( A == 'a') B == 'b',
        if ( A == 'b') x >= 0,
        u + v >= w,
        if ( u <= 0 ) w >= 0
    ))
    b <- blocks(E)
    expect_equal(length(b),2)
    expect_equal(sort(getVars(b[[1]])),c("A","B","x","y","z"))
    expect_equal(sort(getVars(b[[2]])),c("u","v","w"))
})



