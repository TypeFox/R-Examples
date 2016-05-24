
context("Imputation using linear restrictions")


test_that("imputation by pseudoinverse",{
  # example from de Waal et al (2009) pp 304.
  v <- validator(
    x1 + x2 == x3
    , x2 == x4
    , x5 + x6 + x7 == x8
    , x3 + x8 == x9
    , x9 - x10 == x11
  )
  dat <- data.frame(
    x1 = 145
    , x2 = NA
    , x3 = 155
    , x4 = NA
    , x5 = NA
    , x6 = NA
    , x7 = NA
    , x8 = 86
    , x9 = NA
    , x10 = 217
    , x11 = NA
  )
  
  d2 <- data.frame(
    x1 = 145
    , x2 = 10
    , x3 = 155
    , x4 = 10
    , x5 = NA_real_
    , x6 = NA_real_
    , x7 = NA_real_
    , x8 = 86
    , x9 = 241
    , x10 = 217
    , x11 = 24
  )
  lc <- v$linear_coefficients()
  expect_equivalent(
    pivimpute(A=lc$A, b= lc$b,ops = lc$operators, x=t(dat), eps=1e-8)
  , t(d2))
  
})

test_that("imputation with zeros",{
  # example from de Waal et al (2009) pp 307.
  v <- validator(
    x1 + x2 == x3
    , x2 == x4
    , x5 + x6 + x7 == x8
    , x3 + x8 == x9
    , x9 - x10 == x11
    , x6 > 0, x7 > 0
  )
  
  dat <- data.frame(
    x1 = 145
    , x2 = 10
    , x3 = 155
    , x4 = 10
    , x5 = 86
    , x6 = NA_real_
    , x7 = NA_real_
    , x8 = 86
    , x9 = 241
    , x10 = 217
    , x11 = 24
  )
  d2 <- data.frame(
    x1 = 145
    , x2 = 10
    , x3 = 155
    , x4 = 10
    , x5 = 86
    , x6 = 0
    , x7 = 0
    , x8 = 86
    , x9 = 241
    , x10 = 217
    , x11 = 24
  )
  lc <- v$linear_coefficients()
  expect_equivalent(
    zeroimpute(A=lc$A, b=lc$b, ops=lc$operators, x=t(dat) ,eps=1e-8)
    , t(d2)
  )
  
})



test_that("imputation of zeros",{
  v <- validator(
    x1 + x2 + x3 == x4
    , x2 > 0, x3 > 0
  )
  X <- matrix(c(
    1,NA,NA,1,
    NA,1,1,2
  ),ncol=2)
  attr(X,"changed") <- FALSE
  out <- matrix(c(1,0,0,1, NA,1,1,2),ncol=2)
  attr(out,"changed") <- TRUE
  lc <- v$linear_coefficients()
  expect_equal(
    zeroimpute(A=lc$A, b=lc$b, ops=lc$operators, x=X)
    ,  out
  )
})


test_that("imputation of implied values",{
  
  
  v <- validator(x + 2*y == 3,   x + y + z == 7)
  L <- v$linear_coefficients()
  # takes two iterations to impute
  x_ <- c(1,NA,NA)
  expect_equal(
    impute_implied_x(L$A, L$b, L$operators, x_)
    ,c(1,1,5)
  )
  # inequalities (single iteration)
  v <- validator( x <= 0, x >=0, y <=1, y >= 1)
  L <- v$linear_coefficients()
  x_ <- c(NA,NA)
  expect_equal(
      impute_implied_x(L$A, L$b, L$operators, x_)
      , c(0,1)
  )
})

test_that("impute_lr errors when it should",{
  v <- validator(x == y)
  expect_error(impute_lr(data.frame(x="a",y=10),v), regexp="Linear restrictions on nonnumeric data")
  
})


test_that("impute_lr",{
  # example from DCAR 
  v <- validate::validator( 
    x1 + x2      == x3 
    , x4 + x5 + x6 == x1 
    , x7 + x8      == x2 
    , x9 + x10     == x3) 
  dat <- data.frame( 
    x1 = 100, x2=NA_real_, x3=NA_real_, x4 = 15, x5 = NA_real_ 
    , x6 = NA_real_, x7 = 25, x8 = 35, x9 = NA_real_, x10 = 5) 
  
  expect_equal(
    impute_lr(dat,v) 
  , data.frame( x1 = 100, x2=60, x3=160, x4 = 15, x5 = NA_real_ 
              , x6 = NA_real_, x7 = 25, x8 = 35, x9 = 155, x10 = 5)
  ) 
})

test_that("imputation by range determination",{
  # y == 2
  # y + z == 3  (so z=1)
  # x + y <= 0  
  A <- matrix(c(0,0,1, 1,1,1,0,1,0),nrow=3)
  b <- c(2,3,0)
  impute_range_x(x=c(NA,NA,NA), A=A, b=b, neq=2, nleq=1, eps=1e-8)
  lintools::ranges(A,b,neq=2,nleq=1)
})




