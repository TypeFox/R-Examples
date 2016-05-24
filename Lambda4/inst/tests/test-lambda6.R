context("lambda6 Coefficient")
test_that('lambda6 and population covariance matrices',{
  
  cong.1f<-as.numeric(round(lambda6(cong1f)[[1]][[1]],6))
  expect_that(cong.1f, equals(.836908))
  
  par.1f<-as.numeric(round(lambda6(par1f)[[1]][[1]],6))
  expect_that(par.1f, equals(.875))
  
  tau.1f<-as.numeric(round(lambda6(tau1f)[[1]][[1]],6))
  expect_that(tau.1f, equals(.815579))
  
  
  cong.3f<-as.numeric(round(lambda6(cong3f)[[1]][[1]],6))
  expect_that(cong.3f, equals(.785044))
  
  par.3f<-as.numeric(round(lambda6(par3f)[[1]][[1]],6))
  expect_that(par.3f, equals(.832155))
  
  tau.3f<-as.numeric(round(lambda6(tau3f)[[1]][[1]],6))
  expect_that(tau.3f, equals(.760748))
  
  
  cong.5f<-as.numeric(round(lambda6(cong5f)[[1]][[1]],6))
  expect_that(cong.5f, equals(.836421))
  
  par.5f<-as.numeric(round(lambda6(par5f)[[1]][[1]],6))
  expect_that(par.5f, equals(.87367))
  
  tau.5f<-as.numeric(round(lambda6(tau5f)[[1]][[1]],6))
  expect_that(tau.5f, equals(.816616))
  
})