context("variable elimination")




test_that("Gaussian elimination",{
  A <- matrix(c(
    1,1,-1
    ,0,1,-1),byrow=TRUE,nrow=2
  )
  b <- c(1,0)
  neq=2
  
  L <- eliminate(A=A,b=b,neq=neq,variable=1)
  expect_equivalent(L$A,matrix(c(0,1,-1)))
  expect_equivalent(L$b,0)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$H,NULL)
  expect_equal(L$h,0)
})

test_that("Eliminate named variable",{
  A <- matrix(c(
     1,1,-1
    ,0,1,-1),byrow=TRUE,nrow=2
  )
  colnames(A) <- paste0("x",1:3)
  b <- c(1,0)
  neq=2
  
  L <- eliminate(A=A,b=b,neq=neq,variable='x1')
  expect_equivalent(L$A,matrix(c(0,1,-1)))
  expect_equivalent(L$b,0)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$H,NULL)
  expect_equal(L$h,0)
})



test_that("Fourier-Motzkin elimination",{
  
  A <- matrix(c(
    4, -5, -3,  1,
   -1,  1, -1,  0,
    1,  1,  2,  0,
   -1,  0,  0,  0,
    0, -1,  0,  0,
    0,  0, -1,  0),byrow=TRUE,nrow=6) 
  
  b <- c(0,2,3,0,0,0)
  L <- eliminate(A=A,b=b,neq=0,nleq=6,variable=1)
 
  expect_equivalent(L$A,
    matrix(c(
      0, -0.25, -1.75, 0.25,
      0, -1.25, -0.75, 0.25,
      0,  2.00,  1.00, 0.00,
      0,  1.00,  2.00, 0.00,
      0, -1.00,  0.00, 0.00,
      0,  0.00, -1.00, 0.00), byrow=TRUE,nrow=6)
  )
  expect_equivalent(L$b,c(2,0,5,3,0,0))
  expect_equivalent(L$H,
    matrix(c(
      TRUE,  TRUE, FALSE, FALSE, FALSE, FALSE,
      TRUE, FALSE, FALSE,  TRUE, FALSE, FALSE,
      FALSE,  TRUE,  TRUE, FALSE, FALSE, FALSE,
      FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE,  TRUE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE,  TRUE), byrow=TRUE,nrow=6)
  )
  expect_equal(L$h,1)
  expect_equal(L$nleq,6)
  
  # x + y - z <= 0
  #   - y     < 1
  A <- matrix(c(
    1,1,-1,
    0,-1,0
  ),byrow=TRUE,nrow=2)
  b <- c(0,1)
  
  # x - z < 1 ( "<" dominates "<=" )
  L <- eliminate(A,b,neq=0,nleq=1,variable=2)
  expect_equivalent(L$A,matrix(c(1,0,-1),nrow=1))
  expect_equivalent(L$b,1)
  expect_equal(L$neq,0)
  expect_equal(L$nleq,0)
  expect_equal(L$H,matrix(c(TRUE,TRUE),nrow=1))
  expect_equal(L$h,1)
  
})


test_that("elimination with equality and inequalities",{
  
  A <- matrix(c(
     1,   1,  -1,
    -1,   0,   0,
     0,  -1,   0,
     0,   0,  -1),nrow=4,byrow=TRUE)
  
  b <- rep(0,4)
  L <- eliminate(A,b,neq=1,variable=1)
  expect_equivalent(L$A,matrix(c(0,1,-1,0,-1,0,0,0,-1),nrow=3,byrow=TRUE))
  expect_equivalent(L$b,rep(0,3))
  expect_equal(L$neq,0)
  expect_equivalent(L$H,
    matrix(c(
       TRUE,  TRUE, FALSE, FALSE,
      FALSE, FALSE,  TRUE, FALSE,
      FALSE, FALSE, FALSE,  TRUE
    ),nrow=3,byrow=TRUE)
  )
  expect_equal(L$h,1)
  
})



test_that("bugfixes",{
  A <- matrix(c(1,1,-1,0,1,0),nrow=2)
  b <- c(0,0)
  expect_equal(eliminate(A=A,b = b,neq=2,variable = 1)$neq, 1)
  
  A <- matrix(c(
    0,1,0,
    0,1,1,
    1,1,0),byrow=TRUE,nrow=3
  )
  b <- c(2,3,0)
  expect_equal(eliminate(A=A,b=b,neq=2,nleq=1,variable=1)$nleq,0)
})


