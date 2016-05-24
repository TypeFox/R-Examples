library(EcoSimR)
context("Size ratio null model tests")
test_that("size_uniform algorithm works:",{
  expect_true(is.vector(size_uniform()))
  expect_true(is.vector(size_uniform(dataRodents$Sonoran)))
  
  }
)

test_that("size_uniform_user algorithm works:",{
  expect_true(is.vector(size_uniform_user()))
  expect_true(is.vector(size_uniform_user(dataRodents$Sonoran)))
  expect_true(is.vector(do.call(size_uniform_user,list(speciesData = dataRodents$Sonoran,userLow=4,userHigh=20))))
  expect_true(min(do.call(size_uniform_user,list(speciesData = dataRodents$Sonoran,userLow=4,userHigh=20))) > 4)
  expect_true(max(do.call(size_uniform_user,list(speciesData = dataRodents$Sonoran,userLow=4,userHigh=20))) < 20)
  }
)

test_that("size_source_pool algorithm works:",{
  expect_true(is.vector(size_source_pool()))
  expect_true(is.vector(size_source_pool(dataRodents$Sonoran)))
  expect_true(is.vector(do.call(size_source_pool,list(speciesData = dataRodents$Sonoran,
       sourcePool = runif(1000,min(dataRodents$Sonoran),max(dataRodents$Sonoran)),
       speciesProbs = rbeta(1000,1,1) ))))
}
)

test_that("size_size_gamma algorithm works:",{
  expect_true(is.vector(size_gamma()))
  expect_true(is.vector(size_gamma(dataRodents$Sonoran))) 
}
)


test_that("min_diff metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(min_diff()))
  ### Test that minimum difference is accurate
  expect_equal(min_diff(c(1,1.2,-3,3)),.2) 
})

test_that("min_diff metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(min_diff()))
  ### Test that minimum difference is accurate
  expect_equal(min_diff(c(1,1.2,-3,3)),.2) 
  expect_equal(min_diff(c(1,1.2,1.1,-3,3)),.1) 
  
  
  })

test_that("min_ratio metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(min_ratio()))
  ### Test that minimum ratio is accurate
  expect_equal(min_ratio(c(1,2,3,4,5,6)),(6/5)) 
  
})


test_that("var_diff metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(var_diff()))
})


test_that("var_ratio metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(var_ratio()))
  })

test_that("size_null_model works with all combinations of metrics and algorithms",{
  ### Test that proper object is returned
  expect_is(size_null_model(dataRodents,metric ="min_diff" ,algo = "size_uniform",nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="min_ratio" ,algo = "size_uniform",nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_diff" ,algo = "size_uniform",nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_ratio" ,algo = "size_uniform",nRep=10),"sizenullmod")
  
  expect_is(size_null_model(dataRodents,metric ="min_diff" ,algo = "size_uniform_user",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="min_ratio" ,algo = "size_uniform_user",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_diff" ,algo = "size_uniform_user",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_ratio" ,algo = "size_uniform_user",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  
  
  expect_is(size_null_model(dataRodents,metric ="min_diff" ,algo = "size_source_pool",algoOpts = list(sourcePool = runif(1000,min(dataRodents$Sonoran),max(dataRodents$Sonoran)),
                                                                                             speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="min_ratio" ,algo = "size_source_pool",algoOpts = list(sourcePool = runif(1000,min(dataRodents$Sonoran),max(dataRodents$Sonoran)),
                                                                                              speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_diff" ,algo = "size_source_pool",algoOpts = list(sourcePool = runif(1000,min(dataRodents$Sonoran),max(dataRodents$Sonoran)),
                                                                                                      speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_ratio" ,algo = "size_source_pool",algoOpts = list(sourcePool = runif(1000,min(dataRodents$Sonoran),max(dataRodents$Sonoran)),
                                                                                                       speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
   
  expect_is(size_null_model(dataRodents,metric ="min_diff" ,algo = "size_gamma",nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="min_ratio" ,algo = "size_gamma",nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_diff" ,algo = "size_gamma",nRep=10),"sizenullmod")
  expect_is(size_null_model(dataRodents,metric ="var_ratio" ,algo = "size_gamma",nRep=10),"sizenullmod")
  
  smod <- size_null_model(dataRodents,metric ="var_diff" ,algo = "size_gamma",nRep=100)
  
  expect_output(summary(smod),"Metric:  var_diff")
  expect_true(is.null(plot(smod,type="hist")))
  expect_true(is.list(plot(smod,type="size")))
  
  
})

test_that("all text data frames are handled proprely",{
  expect_is(size_null_model(dataRodents,metric ="var_ratio" ,algo = "size_gamma",nRep=10),"sizenullmod")
  
  
  
  
})


