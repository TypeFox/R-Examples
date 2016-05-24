context("Tests bayesian regression code (ictregBayes)")
rm(list=ls())

set.seed(1)

library(MASS)

data(multi)

test_that("ictregBayes works", {
  
  skip_on_cran()
  
  ## Multiple chain MCMC list experiment regression
  ## starts with overdispersed MLE starting values
  
  ## Standard single sensitive-item design
  
  ## Control item parameters fully constrained
  
  mle.estimates <- ictreg(y ~ male + college + age + south, data = race)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
                   Sigma = vcov(mle.estimates) * 9)
  
  bayesDraws.1 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[1, 1:5], psi.start = draws[1, 6:10], burnin = 100,
                              n.draws = 500, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
                              constrained.single = "full")
  
  bayesDraws.2 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[2, 1:5], psi.start = draws[2, 6:10], burnin = 100,
                              n.draws = 500, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
                              constrained.single = "full")
  
  bayesDraws.3 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[3, 1:5], psi.start = draws[3, 6:10], burnin = 100,  
                              n.draws = 500, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
                              constrained.single = "full")
  
  bayesSingleConstrained <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)
  
  summary(bayesSingleConstrained)
  
  ## Control item parameters unconstrained
  
  mle.estimates <- ictreg(y ~ male + college + age + south, data = race,
                          constrained = FALSE)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
                   Sigma = vcov(mle.estimates) * 9)
  
  bayesDraws.1 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[1, 1:5], psi.start = list(psi0 = draws[1, 6:10], 
                                                                            psi1 = draws[1, 11:15]), burnin = 100, n.draws = 500, 
                              delta.tune = diag(.002, 5), 
                              psi.tune = list(psi0 = diag(.0017, 5), psi1 = diag(.00005, 5)), 
                              constrained.single = "none")
  
  bayesDraws.2 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[2, 1:5], psi.start = list(psi0 = draws[2, 6:10], 
                                                                            psi1 = draws[2, 11:15]), burnin = 100, n.draws = 500, 
                              delta.tune = diag(.002, 5), 
                              psi.tune = list(psi0 = diag(.0017, 5), psi1 = diag(.00005, 5)),
                              constrained.single = "none")
  
  bayesDraws.3 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[3, 1:5], psi.start = list(psi0 = draws[3, 6:10], 
                                                                            psi1 = draws[3, 11:15]), burnin = 100, n.draws = 500, 
                              delta.tune = diag(.002, 5), 
                              psi.tune = list(psi0 = diag(.0017, 5), psi1 = diag(.00005, 5)),
                              constrained.single = "none")
  
  bayesSingleUnconstrained <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)
  
  summary(bayesSingleUnconstrained)
  
  ## Control item parameters constrained except intercept
  
  mle.estimates <- ictreg(y ~ male + college + age + south, data = race,
                          constrained = TRUE)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
                   Sigma = vcov(mle.estimates) * 9)
  
  bayesDraws.1 <-  ictregBayes(y ~ male + college + age + south, data = race, 
                               delta.start = draws[1, 1:5], psi.start = c(draws[1, 6:10],0),
                               burnin = 100, n.draws = 500, delta.tune = diag(.002, 5),
                               psi.tune = diag(.0004, 6), constrained.single = "intercept")
  
  bayesDraws.2 <-  ictregBayes(y ~ male + college + age + south, data = race, 
                               delta.start = draws[2, 1:5], psi.start = c(draws[2, 6:10],0),
                               burnin = 100, n.draws = 500, delta.tune = diag(.002, 5),
                               psi.tune = diag(.0004, 6), constrained.single = "intercept")
  
  bayesDraws.3 <-  ictregBayes(y ~ male + college + age + south, data = race, 
                               delta.start = draws[3, 1:5], psi.start = c(draws[3, 6:10],0),
                               burnin = 100, n.draws = 500, delta.tune = diag(.002, 5),
                               psi.tune = diag(.0004, 6), constrained.single = "intercept")
  
  bayesSingleInterceptOnly <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)
  
  summary(bayesSingleInterceptOnly)
  
  ##stop("This gives off a warning it shouldn't.")
  
  ## Multiple sensitive item design
  
  ## Constrained (estimated control item count not included in sensitive fit)
  
  mle.estimates.multi <- ictreg(y ~ male + college + age + south, data = multi,
                                constrained = TRUE)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
                   Sigma = vcov(mle.estimates.multi) * 9)
  
  bayesMultiDraws.1 <- ictregBayes(y ~ male + college + age + south, 
                                   data = multi, delta.start = list(draws[1, 6:10], draws[1, 11:15]), 
                                   psi.start = draws[1, 1:5], burnin = 100, n.draws = 500, 
                                   delta.tune = diag(.002, 5), psi.tune = diag(.001, 5), 
                                   constrained.multi = TRUE)
  
  bayesMultiDraws.2 <- ictregBayes(y ~ male + college + age + south, 
                                   data = multi, delta.start = list(draws[2, 6:10], draws[2, 11:15]), 
                                   psi.start = draws[2, 1:5], burnin = 100, n.draws = 500, 
                                   delta.tune = diag(.002, 5), psi.tune = diag(.001, 5), 
                                   constrained.multi = TRUE)
  
  bayesMultiDraws.3 <- ictregBayes(y ~ male + college + age + south, 
                                   data = multi, delta.start = list(draws[3, 6:10], draws[3, 11:15]), 
                                   psi.start = draws[3, 1:5], burnin = 100, n.draws = 500, 
                                   delta.tune = diag(.002, 5), psi.tune = diag(.001, 5), 
                                   constrained.multi = TRUE)
  
  bayesMultiConstrained <- as.list(bayesMultiDraws.1, bayesMultiDraws.2, 
                                   bayesMultiDraws.3)
  
  summary(bayesMultiConstrained)
  
  ## Unconstrained (estimated control item count is included in sensitive fit)
  
  mle.estimates.multi <- ictreg(y ~ male + college + age + south, data = multi,
                                constrained = FALSE)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
                   Sigma = vcov(mle.estimates.multi) * 9)
  
  bayesMultiDraws.1 <- ictregBayes(y ~ male + college + age + south, 
                                   data = multi, delta.start = list(draws[1, 6:10], draws[1, 11:15]), 
                                   psi.start = draws[1, 1:5], burnin = 100, n.draws = 500, 
                                   delta.tune = diag(.0085, 6), psi.tune = diag(.00025, 5), 
                                   constrained.multi = FALSE)
  
  bayesMultiDraws.2 <- ictregBayes(y ~ male + college + age + south, 
                                   data = multi, delta.start = list(draws[2, 6:10], draws[2, 11:15]), 
                                   psi.start = draws[2, 1:5], burnin = 100, n.draws = 500, 
                                   delta.tune = diag(.0085, 6), psi.tune = diag(.00025, 5), 
                                   constrained.multi = FALSE)
  
  bayesMultiDraws.3 <- ictregBayes(y ~ male + college + age + south, 
                                   data = multi, delta.start = list(draws[3, 6:10], draws[3, 11:15]), 
                                   psi.start = draws[3, 1:5], burnin = 100, n.draws = 500, 
                                   delta.tune = diag(.0085, 6), psi.tune = diag(.00025, 5), 
                                   constrained.multi = FALSE)
  
  bayesMultiUnconstrained <- as.list(bayesMultiDraws.1, bayesMultiDraws.2, 
                                     bayesMultiDraws.3)
  
  summary(bayesMultiUnconstrained)
  
  ##stop("this part gives a warning it shouldn't")
  
  ## Mixed effects models
  
  ## Varying intercepts
  
  mle.estimates <- ictreg(y ~ male + college + age + south, data = race)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
                   Sigma = vcov(mle.estimates) * 9)
  
  bayesDraws.1 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[1, 1:5], psi.start = draws[1, 6:10], burnin = 100, 
                              n.draws = 500, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
                              constrained.single = "full", group.mixed = "state", formula.mixed = ~ 1)
  
  bayesDraws.2 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[2, 1:5], psi.start = draws[2, 6:10], burnin = 100, 
                              n.draws = 500, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
                              constrained.single = "full", group.mixed = "state", formula.mixed = ~ 1)
  
  bayesDraws.3 <- ictregBayes(y ~ male + college + age + south, data = race, 
                              delta.start = draws[3, 1:5], psi.start = draws[3, 6:10], burnin = 100, 
                              n.draws = 500, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
                              constrained.single = "full", group.mixed = "state", formula.mixed = ~ 1)
  
  bayesMixed <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)
  
  summary(bayesMixed)
  
  
})

test_that("predict for bayes works", {
  
  skip_on_cran()
  
  bayes.fit <- ictregBayes(y ~ age + college + male + south, data = multi, 
                           treat = "treat", 
                           delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5),
                           n.draws = 500, burnin = 100, thin = 0)
  
  bayes.predict <- predict(bayes.fit, interval = "confidence", se.fit = TRUE)
  
  
})
