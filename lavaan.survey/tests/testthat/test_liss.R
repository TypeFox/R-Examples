
library(lavaan.survey)

context("LISS example")


test_that("estimate matches", {
  skip_on_cran()
  
  data(liss)
  
  # Fit the model using listwise deletion
  fit.liss <- lavaan("
     cs08 =~ 1 * cs08a247
     cs09 =~ 1 * cs09b247
     cs10 =~ 1 * cs10c247
     cs11 =~ 1 * cs11d247
   
     cs09 ~ cs08
     cs10 ~ cs09
     cs11 ~ cs10
   
     cs08a247 ~~ vare * cs08a247
     cs09b247 ~~ vare * cs09b247
     cs10c247 ~~ vare * cs10c247
     cs11d247 ~~ vare * cs11d247
   
     cs08 ~~ vart08 * cs08
   
     reliab.ratio := vart08 / (vart08 + vare)
  "
  , auto.var = TRUE, meanstructure = TRUE, 
    int.ov.free = TRUE, data = liss)
    
  # Fit the model accounting for nesting of respondents within households
  des.liss <- svydesign(ids = ~nohouse_encr, prob = ~1, data = liss)
  fit.liss.surv <- lavaan.survey(fit.liss, des.liss)
  
  # Complex survey inference on the reliability of interest:
  
  
  ## To deal with missing data (including attrition), multiple imputation can be used.
  ## For example using the mice library (although any MI software is suitable)
  
  ## Uncomment below to run this time-intensive analysis
  ## NOT RUN:
  
  #set.seed(20140221)
  library("mice") 
  #liss.imp <- mice(liss, m = 100, method = "norm", maxit = 100)
  
  ## Turn the mice object into a list() of imputed datasets
  #liss.implist <- lapply(seq(liss.imp$m), function(im) complete(liss.imp, im))
  
  ## After obtaining the list of imputed datasets, 
  ##  use the mitools package to turn it into an imputation list
  library("mitools")
  #liss.implist <- imputationList(liss.implist)
  #save(liss.implist, file="~/Dropbox/Development/lavaan.survey/lavaan.survey/tests/testthat/liss_implist.rdata")
  
  load("liss_implist.rdata")
  
  ## Give the imputation list as data to a svydesign object
  des.liss.imp <- svydesign(ids = ~nohouse_encr, prob = ~1, data = liss.implist)
  
  ## lavaan.survey can be used as usual, using the 
  ##    svydesign object that has an imputation list as data
  ## Standard errors and chi-square tests will account for both the clustering and the 
  ##   imputation uncertainty applying Rubin's rules. 
  fit.liss.surv.mi <- lavaan.survey(fit.liss, des.liss.imp)


  est_surv <- parameterEstimates(fit.liss.surv)[24, ]
  est_surv_mi <- parameterEstimates(fit.liss.surv.mi)[24, ]
  
  expect_equal(est_surv[,'est'], 0.6218929, tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(est_surv[,'se'], 0.01656459, tolerance = 1e-6, check.attributes = FALSE)
  
  expect_equal(est_surv_mi[,'est'], 0.6122555, tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(est_surv_mi[,'se'], 0.01092366, tolerance = 1e-6, check.attributes = FALSE)
})

