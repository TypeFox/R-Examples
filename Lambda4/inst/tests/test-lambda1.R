context("Lambda1 Coefficient")
test_that('lambda1 and population covariance matrices',{
  
  cong.1f<-as.numeric(round(lambda1(cong1f)[[1]],6))
  expect_that(cong.1f, equals(.744627))
  
  par.1f<-as.numeric(round(lambda1(par1f)[[1]],6))
  expect_that(par.1f, equals(.777778))
  
  tau.1f<-as.numeric(round(lambda1(tau1f)[[1]],6))
  expect_that(tau.1f, equals(.729378))
  
  
  cong.3f<-as.numeric(round(lambda1(cong3f)[[1]],6))
  expect_that(cong.3f, equals(.691979))
  
  par.3f<-as.numeric(round(lambda1(par3f)[[1]],6))
  expect_that(par.3f, equals(.729730))
  
  tau.3f<-as.numeric(round(lambda1(tau3f)[[1]],6))
  expect_that(tau.3f, equals(.675234))
  
  
  cong.5f<-as.numeric(round(lambda1(cong5f)[[1]],6))
  expect_that(cong.5f, equals(.764733))
  
  par.5f<-as.numeric(round(lambda1(par5f)[[1]],6))
  expect_that(par.5f, equals(.795918))
  
  tau.5f<-as.numeric(round(lambda1(tau5f)[[1]],6))
  expect_that(tau.5f, equals(.750200))
  
})