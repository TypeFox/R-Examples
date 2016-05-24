
context("Editset")
test_that("editset parses categorical edits",{
    
    v <- expression(
        A %in% c('a','b'),
        B %in% c('c','d'),
        if ( A == 'a') B == 'c'
    )
    E <- editset(v)
    expect_equal(E$num,editmatrix(expression()))
    expect_equal(E$mixnum,editmatrix(expression()))
    expect_equal(E$mixcat,editarray(v))
})

test_that("editset parses numerical edits",{
    v <- expression(x + y == z, 2*x -u == v)
    E <- editset(v)
    expect_equal(E$num,editmatrix(v))
    expect_equal(E$mixnum,editmatrix(expression()))
    expect_equal(E$mixcat,editarray(expression()))
})

test_that("editset parses conditional numeric edits",{
    # test 1: inequalities
    v <- expression( if ( x > 0 ) y > 0 )
    E <- editset(v)
    expect_equal(E$num, editmatrix(expression()))
    expect_equivalent(E$mixnum, editmatrix(expression(x>0,y<=0)))
    expect_equivalent(getArr(E$mixcat),array(c(F,T,F,T),dim=c(1,4)))

    # test 2: with equality in if-statement

    v <- expression( if ( x >= 0 ) y >= 0)
    E <- editset(v)
    expect_equal(E$num, editmatrix(expression()))
    expect_equivalent(E$mixnum,editmatrix(expression(x>=0,y<0)))
    expect_equivalent(getArr(E$mixcat), array(c(F,T,F,T),dim=c(1,4)))

})

test_that("editset parses conditional categorical/numerical edits",{
    # test 1: numerical statement in 'then' clause
    v <- expression(
        A %in% letters[1:2],
        B %in% letters[3:4],
        if ( A == 'a' ) x > 0
    )
    E <- editset(v)
    expect_equal(E$num, editmatrix(expression()))
    expect_equivalent(E$mixnum, editmatrix(expression(x<=0)))
    expect_equal(dim(E$mixcat),c(1,6))
    expect_equivalent(getArr(E$mixcat),array(c(T,F,T,T,F,T),dim=c(1,6)))
    
    # test 2: numerical statement in 'then' clause
    v <- expression(
        A %in% letters[1:2],
        B %in% letters[3:4],
        if ( x > 0 ) A == 'a'
    )
    E <- editset(v)    
    expect_equal(E$num, editmatrix(expression()))
    expect_equivalent(E$mixnum, editmatrix(expression(x>0)))
    expect_equivalent(getArr(E$mixcat), array(c(F,T,T,T,F,T),dim=c(1,6)))

  # throws exception in editrules 2.8.0 (thanks to Alois Haslinger)
  editset(expression(
      x %in% letters[1:2]
    , y %in% letters[3:5]
    , z > 0
    , if(x =='a' && y == 'c') z < 7
  ))

})



## various editset functionalities
test_that("contains finds the right variables in an editset",{
    E <- editset(expression(
        x + y == z,
        if ( s + t >= 6 ) x < 10,
        A %in% letters[1:2],
        if ( A == 'a' )  x > 3
    ))
    expect_equivalent(
        contains(E),
        matrix(c(
            T,T,T,F,F,F,
            T,F,F,T,T,F,
            T,F,F,F,F,T),
            byrow=TRUE,
            nrow=3)
    )

})

test_that("as.editset for cateditmatrix works",{
    E <- cateditmatrix(expression(
         gender %in% c("male", "female")
        , if (pregnant) gender == "female"
    ))
  
  as.editset(E)
})


test_that("simple mixed edit without coefficients is recognized",{
    E <- editset(expression(if ( A %in% c('a','b') ) x > 0 ))
    expect_equal(nedits(E),1)
    expect_equal(length(getVars(E,type="dummy")),1)
})


test_that("Mixed parsing edit containing brackets works",{
  E <- editset("if  ((x  >  0)  &&  (y < 0))  z < y")
})




