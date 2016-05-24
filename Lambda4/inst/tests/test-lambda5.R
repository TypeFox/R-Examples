context("lambda5 Coefficient")
test_that('lambda5 and population covariance matrices',{
  
  cong.1f<-as.numeric(round(lambda5(cong1f)[[1]][[1]],6))
  expect_that(cong.1f, equals(.829847))
  
  par.1f<-as.numeric(round(lambda5(par1f)[[1]][[1]],6))
  expect_that(par.1f, equals(.851271))
  
  tau.1f<-as.numeric(round(lambda5(tau1f)[[1]][[1]],6))
  expect_that(tau.1f, equals(.798297))
  
  
  cong.3f<-as.numeric(round(lambda5(cong3f)[[1]][[1]],6))
  expect_that(cong.3f, equals(.740318))
  
  par.3f<-as.numeric(round(lambda5(par3f)[[1]][[1]],6))
  expect_that(par.3f, equals(.77317))
  
  tau.3f<-as.numeric(round(lambda5(tau3f)[[1]][[1]],6))
  expect_that(tau.3f, equals(.71543))
  
  
  cong.5f<-as.numeric(round(lambda5(cong5f)[[1]][[1]],6))
  expect_that(cong.5f, equals(.789218))
  
  par.5f<-as.numeric(round(lambda5(par5f)[[1]][[1]],6))
  expect_that(par.5f, equals(.81742))
  
  tau.5f<-as.numeric(round(lambda5(tau5f)[[1]][[1]],6))
  expect_that(tau.5f, equals(.770467))
  
})