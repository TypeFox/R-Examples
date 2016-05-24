context("Omega Total Coefficient")
test_that('omega.tot and population covariance matrices',{
  
  cong.1f<-as.numeric(round(omega.tot(cong1f, factors=1)[[1]],6))
  expect_that(cong.1f, equals(.856378))
  
  par.1f<-as.numeric(round(omega.tot(par1f, factors=1)[[1]],6))
  expect_that(par.1f, equals(.888889))
  
  tau.1f<-as.numeric(round(omega.tot(tau1f, factors=1)[[1]],6))
  expect_that(tau.1f, equals(.839789))
  
  
  cong.3f<-as.numeric(round(omega.tot(cong3f, factors=3)[[1]],6))
  expect_that(cong.3f, equals(.826696))
  
  par.3f<-as.numeric(round(omega.tot(par3f, factors=3)[[1]],6))
  expect_that(par.3f, equals(.864865))
  
  tau.3f<-as.numeric(round(omega.tot(tau3f, factors=3)[[1]],6))
  expect_that(tau.3f, equals(.807449))
  
  
  cong.5f<-as.numeric(round(omega.tot(cong5f, factors=5)[[1]],6))
  expect_that(cong.5f, equals(.867708))
  
  par.5f<-as.numeric(round(omega.tot(par5f, factors=5)[[1]],6))
  expect_that(par.5f, equals(.897959))
  
  tau.5f<-as.numeric(round(omega.tot(tau5f, factors=5)[[1]],6))
  expect_that(tau.5f, equals(.852201))
  
})