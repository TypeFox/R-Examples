

context("feasibility")

test_that("is_feasible",{
  # x + y == 1 & x + y == 2
  expect_false( is_feasible(A = matrix(rep(1,4),nrow=2), b = c(1,2), neq = 2) )
  # x >= y & x <= y - 1 
  expect_false( is_feasible(A = matrix(c(-1,1,1,-1),byrow=TRUE,nrow=2), b = c(0,-1), neq = 0) )
  expect_false(is_feasible(matrix(0),b=1,neq=1))
  # x + y == 0
  # x > 0
  # y > 0
  expect_false(is_feasible(A = matrix(c(1,1,-1,0,0,-1),byrow=TRUE,nrow=3),b=c(0,0,0),neq=1,nleq=0))
  # x + y == 0
  # x >= 0
  # y >= 0
  expect_true(is_feasible(A = matrix(c(1,1,-1,0,0,-1),byrow=TRUE,nrow=3),b=c(0,0,0),neq=1,nleq=2))
  
})


test_that("simple contradiction detection",{
  # 0 < -1
  expect_false(has_contradiction(A=matrix(0),b=1e-9,neq=0,nleq=0,eps=1e-8))
  expect_true(has_contradiction(A=matrix(0),b=-1,neq=0,nleq=0,eps=1e-8))
  
  # 0 <= 0 
  expect_false( has_contradiction(A=matrix(0),b=1e-9,neq=0,nleq=1,eps=1e-8) )
  expect_true(has_contradiction(A=matrix(0),b=-2e-8,neq=0,nleq=1,eps=1e-8))
  
  # 0 == 1
  expect_false(has_contradiction(A=matrix(0), b=1e-9,neq=1,nleq=0,eps=1e-8))
  expect_true(has_contradiction(A=matrix(0), b=2e-8,neq=1,nleq=0,eps=1e-8))
  
})


