context("Angoff")
test_that('angoff and population covariance matrices',{
  
  cong.1f<-as.numeric(round(angoff(cong1f)[[1]],6))
  expect_that(cong.1f, equals(.855697))
  
  par.1f<-as.numeric(round(angoff(par1f)[[1]],6))
  expect_that(par.1f, equals(.888889))
  
  tau.1f<-as.numeric(round(angoff(tau1f)[[1]],6))
  expect_that(tau.1f, equals(.833968))
  
  
  cong.3f<-as.numeric(round(angoff(cong3f)[[1]],6))
  expect_that(cong.3f, equals(.825914))
  
  par.3f<-as.numeric(round(angoff(par3f)[[1]],6))
  expect_that(par.3f, equals(.864865))
  
  tau.3f<-as.numeric(round(angoff(tau3f)[[1]],6))
  expect_that(tau.3f, equals(.800821))
  
  
  cong.5f<-as.numeric(round(angoff(cong5f)[[1]],6))
  expect_that(cong.5f, equals(.867068))
  
  par.5f<-as.numeric(round(angoff(par5f)[[1]],6))
  expect_that(par.5f, equals(.897959))
  
  tau.5f<-as.numeric(round(angoff(tau5f)[[1]],6))
  expect_that(tau.5f, equals(.846720))
  
  
})