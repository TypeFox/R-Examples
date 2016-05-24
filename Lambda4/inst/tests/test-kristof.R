context("Kristof Coefficient")
test_that('kristof and population covariance matrices',{
  
  cong.1f<-as.numeric(round(kristof(cong1f)[[1]],6))
  expect_that(cong.1f, equals(.854614))
  
  par.1f<-as.numeric(round(kristof(par1f)[[1]],6))
  expect_that(par.1f, equals(.888889))
  
  tau.1f<-as.numeric(round(kristof(tau1f)[[1]],6))
  expect_that(tau.1f, equals(.833575))
  
  
  cong.3f<-as.numeric(round(kristof(cong3f)[[1]],6))
  expect_that(cong.3f, equals(.801693))
  
  par.3f<-as.numeric(round(kristof(par3f)[[1]],6))
  expect_that(par.3f, equals(.841216))
  
  tau.3f<-as.numeric(round(kristof(tau3f)[[1]],6))
  expect_that(tau.3f, equals(.778395))
  
  
  cong.5f<-as.numeric(round(kristof(cong5f)[[1]],6))
  expect_that(cong.5f, equals(.848932))
  
  par.5f<-as.numeric(round(kristof(par5f)[[1]],6))
  expect_that(par.5f, equals(.880476))
  
  tau.5f<-as.numeric(round(kristof(tau5f)[[1]],6))
  expect_that(tau.5f, equals(.829901))
  
})