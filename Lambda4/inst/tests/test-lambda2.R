context("Lambda2 Coefficient")
test_that('lambda2 and population covariance matrices',{
  
  cong.1f<-as.numeric(round(lambda2(cong1f)[[1]],6))
  expect_that(cong.1f, equals(.853730))
  
  par.1f<-as.numeric(round(lambda2(par1f)[[1]],6))
  expect_that(par.1f, equals(.888889))
  
  tau.1f<-as.numeric(round(lambda2(tau1f)[[1]],6))
  expect_that(tau.1f, equals(.833575))
  
  
  cong.3f<-as.numeric(round(lambda2(cong3f)[[1]],6))
  expect_that(cong.3f, equals(.767955))
  
  par.3f<-as.numeric(round(lambda2(par3f)[[1]],6))
  expect_that(par.3f, equals(.808315))
  
  tau.3f<-as.numeric(round(lambda2(tau3f)[[1]],6))
  expect_that(tau.3f, equals(.747951))
  
  
  cong.5f<-as.numeric(round(lambda2(cong5f)[[1]],6))
  expect_that(cong.5f, equals(.813093))
  
  par.5f<-as.numeric(round(lambda2(par5f)[[1]],6))
  expect_that(par.5f, equals(.845246))
  
  tau.5f<-as.numeric(round(lambda2(tau5f)[[1]],6))
  expect_that(tau.5f, equals(.796694))
  
})