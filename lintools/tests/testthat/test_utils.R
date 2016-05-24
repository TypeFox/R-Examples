
context("utilities")

test_that("all_finite works correctly",{
  expect_true(all_finite(c(1:3)))
  expect_false(all_finite(c(1,2,NA)))
  expect_false(all_finite(c(NA,1,2)))
  expect_false(all_finite(c(NaN,1,2)))
  expect_false(all_finite(c(Inf,1,2)))
  expect_false(all_finite(c(-Inf,1,2)))
  expect_false(all_finite(c(1,2,NA)))
  expect_false(all_finite(c(1,2,NaN)))
  expect_false(all_finite(c(1,2,Inf)))
  expect_false(all_finite(c(1,2,-Inf)))
  
})


test_that("dimension checks",{
  expect_error(check_sys(A=matrix(0),b=c(0,0)))
  expect_error(check_sys(A=matrix(0),b=0,neq=2))
  expect_error(check_sys(A=matrix(0,0,nrow=2),b=0,neq=2))
  expect_error(check_sys(A=matrix(0),b=0,x=c(0,0)))
  expect_error(check_sys(A=matrix("a"),b=0))
  expect_error(check_sys(A=matrix(0),b="a"))
  expect_error(check_sys(A=matrix(0),b=0,x="a"))
  expect_error(check_sys(A=0,b=0))
  expect_error(check_sys(A=matrix(0),b=0,tol=-1))
  expect_error(check_sys(A=matrix(NA_real_),b=0))
  expect_error(check_sys(A=matrix(0),b=NA_real_))
  expect_error(check_sys(A=matrix(0),b=0,x=NA_real_))
  expect_error(check_sys(A=matrix(0),b=0,x=0,tol=NA_real_))
  
})


