context("Tests hierarchical bayesian regression code (ictregBayesHier)")
rm(list=ls())

set.seed(1)

data(race)
data(multi)

library(MASS)

test_that("ictregBayesHier works", {
  
  skip_on_cran()
  
  ## Multiple chain MCMC list experiment regression
  ## starts with overdispersed MLE starting values
  
  ## Multiple item two level hierarchical model - varying intercepts
  
  mle.estimates.multi <- ictreg(y ~ male + college, data = multi,
                                constrained = TRUE)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
                   Sigma = vcov(mle.estimates.multi) * 9)
  
  bayesDraws.1 <- ictregBayesHier(y ~ male + college,
                                  formula.level.2 = ~ 1, 
                                  delta.start.level.1 = list(draws[1, 8:9], draws[1, 2:3], draws[1, 5:6]),
                                  data = multi, treat = "treat",
                                  delta.tune = list(rep(0.005, 2), rep(0.05, 2), rep(0.05, 2)),
                                  alpha.tune = rep(0.001, length(unique(multi$state))),
                                  J = 3, group.level.2 = "state",
                                  n.draws = 500, burnin = 100, thin = 0)
  
  bayesDraws.2 <- ictregBayesHier(y ~ male + college,
                                  formula.level.2 = ~ 1, 
                                  delta.start.level.1 = list(draws[2, 8:9], draws[2, 2:3], draws[2, 5:6]),
                                  data = multi, treat = "treat",
                                  delta.tune = list(rep(0.005, 2), rep(0.05, 2), rep(0.05, 2)),
                                  alpha.tune = rep(0.001, length(unique(multi$state))),
                                  J = 3, group.level.2 = "state",
                                  n.draws = 500, burnin = 100, thin = 0)
  
  bayesDraws.3 <- ictregBayesHier(y ~ male + college,
                                  formula.level.2 = ~ 1, 
                                  delta.start.level.1 = list(draws[3, 8:9], draws[3, 2:3], draws[3, 5:6]),
                                  data = multi, treat = "treat",
                                  delta.tune = list(rep(0.005, 2), rep(0.05, 2), rep(0.05, 2)),
                                  alpha.tune = rep(0.001, length(unique(multi$state))),
                                  J = 3, group.level.2 = "state",
                                  n.draws = 500, burnin = 100, thin = 0)
  
  bayesHierTwoLevel <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)
  
  summary(bayesHierTwoLevel)
  
  ## Multiple item two level hierarchical model - including covariates
  
  mle.estimates.multi <- ictreg(y ~ male + college, data = multi,
                                constrained = TRUE)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
                   Sigma = vcov(mle.estimates.multi) * 9)
  
  bayesDraws.1 <- ictregBayesHier(y ~ male + college,
                                  formula.level.2 = ~ age, 
                                  delta.start.level.1 = list(draws[1, 8:9], draws[1, 2:3], draws[1, 5:6]),
                                  data = multi, treat = "treat",
                                  delta.tune = list(rep(0.005, 2), rep(0.05, 2), rep(0.05, 2)),
                                  alpha.tune = rep(0.001, length(unique(multi$state))),
                                  J = 3, group.level.2 = "state",
                                  n.draws = 500, burnin = 100, thin = 0)
  
  bayesDraws.2 <- ictregBayesHier(y ~ male + college,
                                  formula.level.2 = ~ age, 
                                  delta.start.level.1 = list(draws[2, 8:9], draws[2, 2:3], draws[2, 5:6]),
                                  data = multi, treat = "treat",
                                  delta.tune = list(rep(0.005, 2), rep(0.05, 2), rep(0.05, 2)),
                                  alpha.tune = rep(0.001, length(unique(multi$state))),
                                  J = 3, group.level.2 = "state",
                                  n.draws = 500, burnin = 100, thin = 0)
  
  bayesDraws.3 <- ictregBayesHier(y ~ male + college,
                                  formula.level.2 = ~ age, 
                                  delta.start.level.1 = list(draws[3, 8:9], draws[3, 2:3], draws[3, 5:6]),
                                  data = multi, treat = "treat",
                                  delta.tune = list(rep(0.005, 2), rep(0.05, 2), rep(0.05, 2)),
                                  alpha.tune = rep(0.001, length(unique(multi$state))),
                                  J = 3, group.level.2 = "state",
                                  n.draws = 500, burnin = 100, thin = 0)
  
  bayesHierTwoLevel <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)
  
  summary(bayesHierTwoLevel)
  
})

test_that("predict for ictregBayesHier works", {
  
  skip_on_cran()
    
  mle.estimates.multi <- ictreg(y ~ male + college, data = multi,
                                constrained = TRUE)
  
  draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
                   Sigma = vcov(mle.estimates.multi) * 9)
  
  
  bayes.fit <- ictregBayesHier(y ~ male + college,
                               formula.level.2 = ~ 1, 
                               delta.start.level.1 = list(draws[1, 8:9], draws[1, 2:3], draws[1, 5:6]),
                               data = multi, treat = "treat",
                               delta.tune = list(rep(0.005, 2), rep(0.05, 2), rep(0.05, 2)),
                               alpha.tune = rep(0.001, length(unique(multi$state))),
                               J = 3, group.level.2 = "state",
                               n.draws = 100, burnin = 10, thin = 1)
  
  bayes.predict <- predict(bayes.fit, interval = "confidence", se.fit = TRUE)
  
  
})
