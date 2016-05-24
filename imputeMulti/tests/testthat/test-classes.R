
library(testthat)
library(imputeMulti)

context("Class- mod_imputeMulti")

test_that("new(mod_imputeMulti) works", {
  x <- new("mod_imputeMulti",
           method= "EM",
           mle_iter= 6,
           mle_log_lik= 1000,
           mle_cp= "non.informative",
           mle_x_y= data.frame(x=rnorm(100),y=rnorm(100)))
  
  expect_equal(typeof(x), "S4")
  expect_true(is.mod_imputeMulti(x))
})

test_that("mod_imputeMulti errors", {
  expect_error(new("mod_imputeMulti",
                   method= "abc"))
  
  expect_error(new("mod_imputeMulti",
                   method= "EM",
                   mle_iter= "abc"))
  
  expect_error(new("mod_imputeMulti",
                   method= "EM",
                   mle_iter= -10))
  
  expect_error(new("mod_imputeMulti",
                   method= "EM",
                   mle_iter= 6,
                   mle_log_lik= "abc"))
  
  expect_error(new("mod_imputeMulti",
                   method= "EM",
                   mle_iter= 6,
                   mle_log_lik= 1000,
                   mle_cp= 1234))
  
  expect_error(new("mod_imputeMulti",
                   mle_x_y= list("a", "b", "c")))
  
  expect_error(new("mod_imputeMulti",
                   mle_x_y= matrix(rnorm(100), 10)))
  
  expect_error(new("mod_imputeMulti",
                   mle_x_y= 1:100))
})

###########################################################
context("Class- imputeMulti")

test_that("new(imputeMulti) works", {
  x <- new("mod_imputeMulti",
           method= "EM",
           mle_iter= 6,
           mle_log_lik= 1000,
           mle_cp= "non.informative",
           mle_x_y= data.frame(x=rnorm(100),y=rnorm(100)))
  
  c <- as.call(list("round", quote(A)))
  
  x2 <- new("imputeMulti", x, Gcall= c, 
            data=list(data.frame(matrix(rnorm(100),10)),
                      data.frame(matrix(rnorm(100),10))),
            nmiss= 1000)
  
  expect_equal(typeof(x2), "S4")
  expect_true(is.mod_imputeMulti(x))
  expect_true(is.mod_imputeMulti(x2))
})


test_that("imputeMulti errors", {
  x <- new("mod_imputeMulti",
           method= "EM",
           mle_iter= 6,
           mle_log_lik= 1000,
           mle_cp= "non.informative",
           mle_x_y= data.frame(x=rnorm(100),y=rnorm(100)))
  
  c <- as.call(list("round", quote(A)))
  
  expect_error(new("imputeMulti", x,
                   Gcall= "abc"))
  
  expect_error(new("imputeMulti", x,
                   call= c))
  
  expect_error(new("imputeMulti", x,
                   Gcall= c,
                   data= as.vector(rep(0, 100), mode= "numeric")))
  
  expect_error(new("imputeMulti", x,
                   Gcall= c,
                   data= matrix(rnorm(100), 10)))
  
  expect_error(new("imputeMulti", x,
                   Gcall= c,
                   data= replicate(5, data.frame(x= rnorm(10))),
                   nmiss= "abc"))
})