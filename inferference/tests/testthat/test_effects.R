context("Effect retrieval functions")

test_that("retrieval functions work", {
  testdt <- subset(vaccinesim, group %in% 1:10)
  allos  <- c(.35, .4)
  
  test <- interference(data = testdt,
                       propensity_integrand = 'logit_integrand',
                       formula = y | A | B ~ X1 + (1|group) | group,
                       allocations = allos,
                       model_method = 'glmer',
                       causal_estimation_options = list(variance_estimation = 'robust'),
                       method = 'simple')
  values <- c('estimate', 'std.error', 'conf.low', 'conf.high')
  # Direct effects
  expect_equivalent(as.numeric(as.matrix(direct_effect(test, .4)[values])), 
                    c(0.24673529536089, 0.122468329766916, 
                      0.00670177977095976, 0.48676881095082))
  # Indirect effects
  expect_equivalent(as.numeric(as.matrix(ie(test, .35, .4)[values])), 
                    c(0.0666863333283656, 0.0679513395569612, 
                      -0.0664958449045302, 0.199868511561261))
  
  # Total effects
  expect_equivalent(as.numeric(as.matrix(te(test, .35, .4)[values])), 
                    c(0.313421628689255  , 0.14769062385725, 
                      0.0239533250747932, 0.602889932303717))
  
  # Overall effects
  expect_equivalent(as.numeric(as.matrix(oe(test, .35, .4)[values])), 
                    c(0.094375104083545  , 0.0714379954812212, 
                      -0.0456407941873837, 0.234391002354474))
})

