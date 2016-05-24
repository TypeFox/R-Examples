context("Raju Coefficient")
test_that('raju and population covariance matrices',{
  
  cong.1f<-as.numeric(round(raju(cong1f)[[1]],6))
  expect_that(cong.1f, equals(.849558))
  
  par.1f<-as.numeric(round(raju(par1f)[[1]],6))
  expect_that(par.1f, equals(.888889))
  
  tau.1f<-as.numeric(round(raju(tau1f)[[1]],6))
  expect_that(tau.1f, equals(.833575))
  
  
  cong.3f<-as.numeric(round(raju(cong3f)[[1]],6))
  expect_that(cong.3f, equals(.819762))
  
  par.3f<-as.numeric(round(raju(par3f)[[1]],6))
  expect_that(par.3f, equals(.864865))
  
  tau.3f<-as.numeric(round(raju(tau3f)[[1]],6))
  expect_that(tau.3f, equals(.800278))
  
  
  cong.5f<-as.numeric(round(raju(cong5f)[[1]],6))
  expect_that(cong.5f, equals(.860936))
  
  par.5f<-as.numeric(round(raju(par5f)[[1]],6))
  expect_that(par.5f, equals(.897959))
  
  tau.5f<-as.numeric(round(raju(tau5f)[[1]],6))
  expect_that(tau.5f, equals(.84638))
  
  
})