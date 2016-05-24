context("Test sinc function")
test_that("Test sinc function", {


  x=seq(-4*pi,4*pi,by=pi/10)


  yx21=sinc(x[21])
  yx41= sinc(x[41])
  yx61= sinc(x[61])

  expect_true( abs(yx21) < 1e-10, info="sinc(n*pi) = 0")
  expect_true( abs(yx40) == 1, info="sinc(0) = 1")
  expect_true( abs(yx61) < 1e-10, info="sinc(n*pi) = 0")

} )

context("Test rcpp hello")
test_that("Test rccp hello", {


 z=rcpp_hello()
 x1=as.character(unlist(z[[1]]))
 x2=as.numeric(unlist(z[[2]]))


  expect_true( x1[1] == "foo", info="x1[1] == foo")
  expect_true( x1[2] == "bar", info="x1[2] == bar")
  expect_true( x2[1] == 0 , info="x2[1] == 0")
  expect_true( x2[2] == 1 , info="x2[2] == 1")

} )
