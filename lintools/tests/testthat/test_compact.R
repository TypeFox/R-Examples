

context("compactify linear systems")

test_that("column removal, no x",{
  L <- compact(
    A = matrix(c(1,0),nrow=1)
    ,b = 1
    ,neq=1
    ,nleq=0
  )
  expect_equivalent(L$A,matrix(1))
  expect_equal(L$b,1)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,TRUE))
  
})

test_that("column removal, with x",{
  L <- compact(
    A = matrix(c(1,0),nrow=1)
    , x = c(2,8)
    ,b = 1
    ,neq=1
    ,nleq=0
  )
  expect_equivalent(L$A,matrix(1))
  expect_equal(L$b,1)
  expect_equal(L$x,2)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,TRUE))
})

test_that("row removal",{
  # x + y == 1
  # 0 + 0 == 0
  L <- compact(
    A = matrix(c(1,1,0,0),nrow=2,byrow=TRUE)
    , b = c(1,0)
    , neq=2
  )
  expect_equivalent(L$A,matrix(c(1,1),nrow=1))
  expect_equal(L$b,1)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,FALSE))
})

test_that("row removal, case w/nonzero b",{
  L <- compact(
    A = matrix(c(1,1,0,0),nrow=2,byrow=TRUE)
    , b = c(1,2)
    , neq=2
  )
  expect_equivalent(L$A,matrix(c(1,1,0,0),nrow=2,byrow=TRUE))
  expect_equal(L$b,c(1,2))
  expect_equal(L$neq,2)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,FALSE))
})


test_that("Combine inequalities, simple case",{
  # x <= 0
  # x >= 0
  L <- compact(
    A=matrix(c(1,-1))
    , b=c(0,0)
    , neq=0
    , nleq=2
    )
 expect_equivalent(L$A, matrix(1,nrow=1))
 expect_equal(L$b,0)
 expect_equal(L$neq,1)
 expect_equal(L$nleq,0)
 expect_equal(L$cols_removed,FALSE)
})

test_that("Combine inequalities, simple case with non-zero b",{
  # x <= 3
  # x >= 3
  L <- compact(
    A=matrix(c(1,-1))
    , b=c(3,3)
    , neq=0
    , nleq=2
    )
 expect_equivalent(L$A, matrix(1,nrow=1))
 expect_equal(L$b,3)
 expect_equal(L$neq,1)
 expect_equal(L$nleq,0)
 expect_equal(L$cols_removed,FALSE)
})


test_that("Combine inequalities, simple case that includes equalities",{
  # x + y == 1
  # x <= 3
  # x >= 3
  L <- compact(
    A=matrix(c(1,1,1,0,-1,0),nrow=3,byrow=TRUE)
    , b=c(1,3,3)
    , neq=1
    , nleq=2
    )
 expect_equivalent(L$A, matrix(c(1,1,1,0),nrow=2,byrow=TRUE))
 expect_equal(L$b,c(1,3))
 expect_equal(L$neq,2)
 expect_equal(L$nleq,0)
 expect_equal(L$cols_removed,c(FALSE,FALSE))
})


test_that("combined inequalities and row removal",{
  #  x + y == 1
  #  x + 0 <= 1
  # -x + 0 <= 1
  #  0 + 0 <=0
  L <- compact(
    A = matrix(c(1,1,1,0,-1,0,0,0),nrow=4,byrow=TRUE)
    , b = c(1,1,1,0)
    , neq = 1
    , nleq = 3
  )
  expect_equivalent(L$A, matrix(c(1,1,1,0),nrow=2,byrow=TRUE))
  expect_equal(L$neq,2)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,FALSE))
})

test_that("earlier bugs",{
  
  A <- matrix(c(1,1,-1,0,0,-1),byrow=TRUE,nrow=3)
  b <- rep(0,3)
  expect_equal(compact(A=A,b=b,neq=1,nleq=2,eps=1e-8)$neq, 1)
  
})


