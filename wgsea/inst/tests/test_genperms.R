##  genperms.R
#it should return a matrix with a length the same as the argument
context("genperms.R\n")

y<-rbinom(50,2,0.3)
gp<-genperms(y,5)
test_that("genperms returns a matrix",{
  expect_that(gp,is_a("matrix"))
})

test_that("genperms returns a matrix with correct dimensions",{
  expect_equal(nrow(gp),50)
  expect_equal(ncol(gp),5)
})
