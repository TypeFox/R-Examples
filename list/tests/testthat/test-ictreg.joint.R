context("Tests list-as-an-outcome ictreg.joint")
rm(list=ls())

set.seed(1)
data(mexico)

test_that("ictreg.joint works", {
  
  skip_on_cran()

  loyal <- mexico[mexico$mex.loyal == 1,]
  notloyal <- mexico[mexico$mex.loyal == 0,]
  
  ## Logistic outcome regression
  ## (effect of vote-selling on turnout)
  ## This replicates Table 4 in Imai et al. 2014
  
  loyalreg <- ictreg.joint(formula = mex.y.all ~ mex.male + mex.age + mex.age2 + mex.education +  
                             mex.interest + mex.married +
                             mex.wealth + mex.urban + mex.havepropoganda + mex.concurrent, data = loyal,
                           treat = "mex.t", outcome = "mex.votecard", J = 3, constrained = TRUE,
                           outcome.reg = "logistic", maxIter = 1000)
  summary(loyalreg)
  
  ## Linear outcome regression
  ## (effect of vote-selling on candidate approval)
  ## This replicates Table 5 in Imai et al. 2014
  
  approvalreg <- ictreg.joint(formula = mex.y.all ~ mex.male + mex.age + mex.age2 +
                                mex.education +
                                mex.interest + mex.married +
                                mex.urban + 
                                mex.cleanelections + mex.cleanelectionsmiss +
                                mex.havepropoganda +
                                mex.wealth + mex.northregion +
                                mex.centralregion + mex.metro + mex.pidpriw2 + 
                                mex.pidpanw2 + mex.pidprdw2,
                              data = mexico, treat = "mex.t", outcome = "mex.epnapprove",
                              J = 3, constrained = TRUE,
                              outcome.reg = "linear", maxIter = 1000)
  
  
  summary(approvalreg)
  
  loyalpred <- predict.ictreg.joint(loyalreg, se.fit = TRUE, interval = "confidence", 
                                    level = 0.95, avg = TRUE, 
                                    sensitive.value = "both", 
                                    sensitive.diff = TRUE, return.draws = TRUE,
                                    predict.sensitive = TRUE)
  
  loyalpred$fit
  
  ## View predicted probability of vote selling, in the sample of party supporters.
  ## This replicates the results in the lefthand panel of Figure 2 in Imai et al. 2014
  
  loyalpred$fitsens
  
})

