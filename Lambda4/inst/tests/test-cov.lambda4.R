context("Covariance Lambda4")
test_that('Covariance Lambda4 and population covariance matrices',{
  
  cong.1f<-as.numeric(round(cov.lambda4(cong1f)[[1]][[1]],6))
  expect_that(cong.1f, equals(.854614))
  
  par.1f<-as.numeric(round(cov.lambda4(par1f)[[1]][[1]],6))
  expect_that(par.1f, equals(.888889))
  
  tau.1f<-as.numeric(round(cov.lambda4(tau1f)[[1]][[1]],6))
  expect_that(tau.1f, equals(.833575))
  
  
  cong.3f<-as.numeric(round(cov.lambda4(cong3f)[[1]][[1]],6))
  expect_that(cong.3f, equals(.823117))
  
  par.3f<-as.numeric(round(cov.lambda4(par3f)[[1]][[1]],6))
  expect_that(par.3f, equals(.864865))
  
  tau.3f<-as.numeric(round(cov.lambda4(tau3f)[[1]][[1]],6))
  expect_that(tau.3f, equals(.800278))
  
  
  cong.5f<-as.numeric(round(cov.lambda4(cong5f)[[1]][[1]],6))
  expect_that(cong.5f, equals(.864896))
  
  par.5f<-as.numeric(round(cov.lambda4(par5f)[[1]][[1]],6))
  expect_that(par.5f, equals(.897959))
  
  tau.5f<-as.numeric(round(cov.lambda4(tau5f)[[1]][[1]],6))
  expect_that(tau.5f, equals(.846380))
  
})