

context('contains')
test_that("contains works for editset",{
  E <- editset(expression(
     if ( x > 0 ) y > 0,
     x >= 0,
     x+y == z
  )) 
  expect_equivalent(
     contains(E,'z'),
     matrix(c(FALSE,TRUE,FALSE),nrow=3,ncol=1)
  )
  # test generated from bug report of Jeroen Pannekoek & MvdL
  E <- editset(expression(  
     0 < v100 + v37 + v38 + v39 + v40 + v41 + v42
    ,0 <= v40
    ,if( 0 < v40 ) v50 >= 1
    ,if( 0 < v40 ) v51 >= 1
    ,if( 0 < v50 ) v40 > 0
    ,if( v40 <= 0 ) 0 >= v51
    ,if( v40 <= 0 ) 0 >= v132
    ,if( 0 < v40 ) v132 >= 1
    ,if( v115 < 1 ) 0 >= v40
  ))
  expect_equivalent(
    contains(E,'v40'),
    matrix(rep(TRUE,9),nrow=9,ncol=1)
  )
})

