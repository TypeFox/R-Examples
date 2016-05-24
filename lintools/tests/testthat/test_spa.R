
context("dense project")

test_that("project",{
  # system:
  # x + y == 1
  # x >  0  ==> -x <= 0
  # y >= 0  ==> -y < 0
  A <- matrix(c(
    1, 1,
    -1, 0,
    0,-1
  ), nrow=3,byrow=TRUE
  )
  b <- c(1,0,0)
  
  out <- project(c(0,0),A,b,neq=1)
  expect_equivalent(out$x , c(0.5,0.5))
  expect_equivalent(out$objective,1/sqrt(2))
  expect_error( project(c(0,NA),A,b,neq=1) )
  expect_error( adjust(c(0,0,0),A,b,neq=1) )
  
  # exceeding max iterations must be noted
  # x < -1, x > 0
  A <- matrix(c(1,-1),nrow=2)
  b <- matrix(c(-1,0))
  expect_equal(project(0,A,b,neq=0)$status,3)
})


context("sparse project")

test_that("sparse_project",{
  A <- data.frame(
    row = c(1,1,2,2)
    ,col= c(1,2,1,2)
    ,coef=c(-1,0,0,-1)
  )
  b <- c(0,0)
  expect_equal(sparse_project(c(-1,-1),A=A,b=b,neq=0)$x,c(0,0))
  expect_equal(sparse_project(c(-1,-1),A=A,b=b,neq=0)$objective, sqrt(2))
})

test_that("sparse_constraints object",{
  A <- data.frame(
    row = c(1,1,2,2,3,3)
    ,col=c(1,2,1,2,1,2)
    ,coef=c(1,1,-1,0,0,-1)
  )
  b<-c(1,0,0)
  sc <- sparse_constraints(A,b,neq=1)
  expect_equal(sc$.nvar(), 2)
  expect_equal(sc$.nconstr(), 3)
  x <- c(1,-2)
  expect_equivalent(sc$.multiply(x),c(-1,-1,2))
  expect_equivalent(sc$.diffvec(x), c(-2,-1,2))
  expect_equivalent(sc$.diffmax(x), 2)
  expect_equivalent(sc$.diffsum(x), 4)
  x <- c(0.5,0.5)
  # no adjusting necessary
  expect_equal(sc$project(x,w=c(1,1),eps=0.01,maxiter=100)$x, x,tolerance=0.01)
  # adjusting necessary
  x <- c(0,0)
  expect_equal(sc$project(x,w=c(1,1), eps=0.01, maxiter=100)$x, c(0.5,0.5), tolerance=0.01)
  
  # no-crash test for printing
  capture.output(print(sc))
})





