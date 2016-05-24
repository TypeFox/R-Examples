context("Lambda3 Coefficient")
test_that('lambda3 and population covariance matrices',{
  
  cong.1f<-as.numeric(round(lambda3(cong1f)[[1]][[1]],6))
  expect_that(cong.1f, equals(.851002))
  
  par.1f<-as.numeric(round(lambda3(par1f)[[1]][[1]],6))
  expect_that(par.1f, equals(.888889))
  
  tau.1f<-as.numeric(round(lambda3(tau1f)[[1]][[1]],6))
  expect_that(tau.1f, equals(.833575))
  
  
  cong.3f<-as.numeric(round(lambda3(cong3f)[[1]][[1]],6))
  expect_that(cong.3f, equals(.754886))
  
  par.3f<-as.numeric(round(lambda3(par3f)[[1]][[1]],6))
  expect_that(par.3f, equals(.796069))
  
  tau.3f<-as.numeric(round(lambda3(tau3f)[[1]][[1]],6))
  expect_that(tau.3f, equals(.736619))
  
  
  cong.5f<-as.numeric(round(lambda3(cong5f)[[1]][[1]],6))
  expect_that(cong.5f, equals(.804982))
  
  par.5f<-as.numeric(round(lambda3(par5f)[[1]][[1]],6))
  expect_that(par.5f, equals(.837809))
  
  tau.5f<-as.numeric(round(lambda3(tau5f)[[1]][[1]],6))
  expect_that(tau.5f, equals(.789685))
  
})