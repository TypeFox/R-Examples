context("Tests regression code (ictreg)")
rm(list=ls())

set.seed(1)

data(race)
data(affirm)
data(multi)

race$weight.1 <- 1
race$weight.uniform <- runif(nrow(race))

affirm$weight.1 <- 1
affirm$weight.uniform <- runif(nrow(affirm))

test_that("diff in means works", {
  
  # Calculate list experiment difference in means
  
  diff.in.means.results <- ictreg(y ~ 1, data = race, treat = "treat", J=3, method = "lm")
  
  expect_is(diff.in.means.results, "ictreg")
  
  expect_equivalent(round(coef(diff.in.means.results), 3), c(.068, 2.134))
  expect_equivalent(round(sqrt(diag(vcov(diff.in.means.results))), 3), c(.05, .033))
  
  diff.in.means.results.w1 <- ictreg(y ~ 1, data = race, treat = "treat", J=3, method = "lm", weights = "weight.1")
  
  expect_that(coef(diff.in.means.results), equals(coef(diff.in.means.results.w1)))
  expect_that(vcov(diff.in.means.results), equals(vcov(diff.in.means.results.w1)))
  
  diff.in.means.results.wU <- ictreg(y ~ 1, data = race, treat = "treat", J=3, method = "lm", weights = "weight.uniform")
  
  expect_that(coef(diff.in.means.results), not(equals(coef(diff.in.means.results.wU))))
  expect_that(vcov(diff.in.means.results), not(equals(vcov(diff.in.means.results.wU))))
  
})

test_that("lm works", {
  
  # Fit linear regression
  # Replicates Table 1 Columns 1-2 Imai (2011); note that age is divided by 10
  
  lm.results <- ictreg(y ~ south + age + male + college, data = race, 
                       treat = "treat", J=3, method = "lm")
  
  expect_is(lm.results, "ictreg")
  
  ## note rescaling of age, third covariate
  sensitive.coef <- c(-.434, .202, .007*10, .180, .114) 
  control.coef <- c(2.406, -.18, .002*10, -.202, -.394)
  
  sensitive.se <- c(.16, .118, .003*10, .098, .098)
  control.se <- c(.105, .074, .002*10, .065, .064)
  
  expect_equivalent(round(coef(lm.results), 2), round(c(sensitive.coef, control.coef), 2))
  expect_equivalent(round(sqrt(diag(vcov(lm.results))), 1), round(c(sensitive.se, control.se), 1))
  
  lm.results.w1 <- ictreg(y ~ south + age + male + college, data = race, 
                          treat = "treat", J=3, method = "lm", weights = "weight.1")
  
  expect_that(coef(lm.results), equals(coef(lm.results.w1)))
  expect_that(vcov(lm.results), equals(vcov(lm.results.w1)))
  
  lm.results.wU <- ictreg(y ~ south + age + male + college, data = race, 
                          treat = "treat", J=3, method = "lm", weights = "weight.uniform")
  
  expect_that(coef(lm.results), not(equals(coef(lm.results.wU))))
  expect_that(vcov(lm.results), not(equals(vcov(lm.results.wU))))
  
})

test_that("nls works", {
  
  # Fit two-step non-linear least squares regression
  # Replicates Table 1 Columns 3-4 Imai (2011); note that age is divided by 10
  
  nls.results <- ictreg(y ~ south + age + male + college, data = race, 
                        treat = "treat", J=3, method = "nls")
  
  expect_is(nls.results, "ictreg")
  
  ## note rescaling of age, third covariate
  sensitive.coef <- c(-7.084, 2.490, .026*10, 3.097, .612)
  control.coef <- c(1.388, -.277, .003*10, -.332, -.662)
  
  sensitive.se <- c(3.669, 1.268, .031*10, 2.929, 1.03)
  control.se <- c(.187, .116, .004*10, .107, .113)
  
  expect_equivalent(round(coef(nls.results), 2), round(c(sensitive.coef, control.coef), 2))
  expect_equivalent(round(sqrt(diag(vcov(nls.results))), 0), round(c(sensitive.se, control.se), 0))
  
  nls.results.w1 <- ictreg(y ~ south + age + male + college, data = race, 
                           treat = "treat", J=3, method = "nls", weights = "weight.1")
  
  expect_that(coef(nls.results), equals(coef(nls.results.w1)))
  expect_that(vcov(nls.results), equals(vcov(nls.results.w1)))
  
  nls.results.wU <- ictreg(y ~ south + age + male + college, data = race, 
                           treat = "treat", J=3, method = "nls", weights = "weight.uniform")
  
  expect_that(coef(nls.results), not(equals(coef(nls.results.wU))))
  expect_that(vcov(nls.results), not(equals(vcov(nls.results.wU))))
  
})

test_that("ml constrained works", {
  
  skip_on_cran()
  
  # Fit EM algorithm ML model with constraint
  # Replicates Table 1 Columns 5-6, Imai (2011); note that age is divided by 10
  
  ml.constrained.results <- ictreg(y ~ south + age + male + college, data = race, 
                                   treat = "treat", J=3, method = "ml", 
                                   overdispersed = FALSE, constrained = TRUE)
  
  expect_is(ml.constrained.results, "ictreg")
  
  ## note rescaling of age, third covariate
  sensitive.coef <- c(-5.508, 1.675, .064*10, .846, -.315)
  control.coef <- c(1.191, -.292, .003*10, -.251, -.516)
  
  sensitive.se <- c(1.021, .559, .016*10, .494, .474)
  control.se <- c(.144, .097, .003*10, .082, .084)
  
  expect_equivalent(round(coef(ml.constrained.results), 2), round(c(sensitive.coef, control.coef), 2))
  expect_equivalent(round(sqrt(diag(vcov(ml.constrained.results))), 1), round(c(sensitive.se, control.se), 1))
  
  ml.constrained.results.w1 <- ictreg(y ~ south + age + male + college, data = race, 
                                      treat = "treat", J=3, method = "ml", 
                                      overdispersed = FALSE, constrained = TRUE, weights = "weight.1")
  
  expect_that(coef(ml.constrained.results), equals(coef(ml.constrained.results.w1)))
  expect_that(vcov(ml.constrained.results), equals(vcov(ml.constrained.results.w1)))
  
  ml.constrained.results.wU <- ictreg(y ~ south + age + male + college, data = race, 
                                      treat = "treat", J=3, method = "ml", 
                                      overdispersed = FALSE, constrained = TRUE, weights = "weight.uniform")
  
  expect_that(coef(ml.constrained.results), not(equals(coef(ml.constrained.results.wU))))
  expect_that(vcov(ml.constrained.results), not(equals(vcov(ml.constrained.results.wU))))
  
})

test_that("ml unconstrained works", {
  
  skip_on_cran()
  
  # Fit EM algorithm ML model with no constraint
  # Replicates Table 1 Columns 7-10, Imai (2011); note that age is divided by 10
  
  ml.unconstrained.results <- ictreg(y ~ south + age + male + college, data = race, 
                                     treat = "treat", J=3, method = "ml", 
                                     overdispersed = FALSE, constrained = FALSE)
  
  expect_is(ml.unconstrained.results, "ictreg")
  
  ## note rescaling of age, third covariate
  sensitive.coef <- c(-6.226, 1.379, .065*10, 1.366, -.182)
  control.coef.psi0 <- c(1.156, -.299, .003*10, -.218, -.488)
  control.coef.psi1 <- c(3.781, -.270, -.013*10, -1.689, -.954)
  
  sensitive.se <- c(1.045, .82, .021*10, .612, .569)
  control.se.psi0 <- c(.156, .107, .003*10, .086, .087)
  control.se.psi1 <- c(2.159, .590, .016*10, 1.633, .715)
  
  expect_equivalent(round(coef(ml.unconstrained.results), 0), 
                    round(c(sensitive.coef, control.coef.psi0, control.coef.psi1), 0))
  expect_equivalent(round(sqrt(diag(vcov(ml.unconstrained.results))), 1), 
                    round(c(sensitive.se, control.se.psi0, control.se.psi1), 1))
  
  ml.unconstrained.results.w1 <- ictreg(y ~ south + age + male + college, data = race, 
                                        treat = "treat", J=3, method = "ml", 
                                        overdispersed = FALSE, constrained = FALSE, weights = "weight.1")
  
  expect_that(coef(ml.unconstrained.results), equals(coef(ml.unconstrained.results.w1)))
  expect_that(vcov(ml.unconstrained.results), equals(vcov(ml.unconstrained.results.w1)))
  
  ml.unconstrained.results.wU <- ictreg(y ~ south + age + male + college, data = race, 
                                        treat = "treat", J=3, method = "ml", 
                                        overdispersed = FALSE, constrained = FALSE, weights = "weight.uniform")
  
  expect_that(coef(ml.unconstrained.results), not(equals(coef(ml.unconstrained.results.wU))))
  expect_that(vcov(ml.unconstrained.results), not(equals(vcov(ml.unconstrained.results.wU))))
  
})

test_that("multi works", {
  
  skip_on_cran()
  
  # Fit EM algorithm ML model for multiple sensitive items
  # Replicates Table 3 in Blair and Imai (2010)
  
  multi.results <- ictreg(y ~ male + college + age + south + south:age, treat = "treat", 
                          J = 3, data = multi, method = "ml", 
                          multi.condition = "level")
  
  expect_is(multi.results, "ictreg")
  
  ## actual estimates from paper
  coef.blackfamily <- c(-7.575, 1.200, -.259, .852, 4.751, -.643, .267)
  se.blackfamily <- c(1.539, .569, .496, .220, 1.85, .347, .252)
  
  coef.affirm <- c(-5.27, .538, -.552, .579, 5.66, -.833, .991)
  se.affirm <- c(1.268, .435, .399, .147, 2.429, .418, .264)
  
  coef.control <- c(1.389, -.325, -.533, .006, -.685, .093)
  se.control <- c(.143, .076, .074, .028, .297, .061)
  
  expect_equivalent(round(coef(multi.results), 1), round(c(coef.affirm, coef.blackfamily, coef.control), 1))
  expect_equivalent(round(sqrt(diag(vcov(multi.results))), 1), round(c(se.control, se.affirm, se.blackfamily), 1))
  
})

test_that("ml model with affirm data works", {
  
  skip_on_cran()
  
  # Fit standard design ML model
  # Replicates Table 7 Columns 1-2 in Blair and Imai (2010)
  
  noboundary.results <- ictreg(y ~ age + college + male + south, treat = "treat",
                               J = 3, data = affirm, method = "ml", 
                               overdispersed = FALSE)
  
  expect_is(noboundary.results, "ictreg")
  
  sensitive.coef <- c(-1.299, .295, -.343, .040, 1.177)
  sensitive.se <- c(.556, .101, .336, .346, .480)
  
  expect_equivalent(round(coef(noboundary.results)[1:5], 2), round(c(sensitive.coef), 2))
  expect_equivalent(round(sqrt(diag(vcov(noboundary.results)))[1:5], 2), round(c(sensitive.se), 2))
  
})

test_that("ceiling model works", {
  
  skip_on_cran()
  
  # Fit standard design ML model with ceiling effects alone
  # Replicates Table 7 Columns 3-4 in Blair and Imai (2010)
  
  ceiling.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
                            J = 3, data = affirm, method = "ml", fit.start = "nls",
                            ceiling = TRUE, ceiling.fit = "bayesglm",
                            ceiling.formula = ~ age + college + male + south)
  
  expect_is(ceiling.results, "ictreg")
  
  sensitive.coef <- c(-1.291, .294, -.345, .038, 1.175)
  sensitive.se <- c(.558, .101, .336, .346, .480)
  
  expect_equivalent(round(coef(ceiling.results)[6:10], 1), round(c(sensitive.coef), 1))
  expect_equivalent(round(sqrt(diag(vcov(ceiling.results)))[6:10], 2), round(c(sensitive.se), 2))
  
})

test_that("floor model works", {
  
  skip_on_cran()
  
  # Fit standard design ML model with floor effects alone
  # Replicates Table 7 Columns 5-6 in Blair and Imai (2010)
  
  floor.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
                          J = 3, data = affirm, method = "ml", fit.start = "glm", 
                          floor = TRUE, floor.fit = "bayesglm",
                          floor.formula = ~ age + college + male + south)
  
  expect_is(floor.results, "ictreg")
  
  sensitive.coef <- c(-1.251, .314, -.605, -.088, .682)
  sensitive.se <- c(.501, .092, .298, .3, .335)
  
  expect_equivalent(round(coef(floor.results)[6:10], 1), round(c(sensitive.coef), 1))
  expect_equivalent(round(sqrt(diag(vcov(floor.results)))[6:10], 2), round(c(sensitive.se), 2))
  
})

test_that("ceiling and floor model works", {
  
  skip_on_cran()
  
  # Fit standard design ML model with floor and ceiling effects
  # Replicates Table 7 Columns 7-8 in Blair and Imai (2010)
  
  both.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
                         J = 3, data = affirm, method = "ml", 
                         floor = TRUE, ceiling = TRUE, 
                         floor.fit = "bayesglm", ceiling.fit = "bayesglm",
                         floor.formula = ~ age + college + male + south,
                         ceiling.formula = ~ age + college + male + south)
  
  expect_is(both.results, "ictreg")
  
  sensitive.coef <- c(-1.245, .313, -.606, -.088, .681)
  sensitive.se <- c(.502, .092, .298, .3, .335)
  
  expect_equivalent(round(coef(both.results)[6:10], 1), round(c(sensitive.coef), 1))
  expect_equivalent(round(sqrt(diag(vcov(both.results)))[6:10], 2), round(c(sensitive.se), 2))
  
})

test_that("plot.predict works for ictreg", {
  
  skip_on_cran()
  
  race.south <- race.nonsouth <- race
  race.south[, "south"] <- 1
  race.nonsouth[, "south"] <- 0
  
  # Fit EM algorithm ML model with constraint
  ml.constrained.results <- ictreg(y ~ south + age + male + college, 
                                   data = race, treat = "treat", J=3, method = "ml", 
                                   overdispersed = FALSE, constrained = TRUE)
  
  # Calculate average predictions for respondents in the South 
  # and the the North of the US for the MLE model, replicating the 
  # estimates presented in Figure 1, Imai (2011)
  avg.pred.south.mle <- predict(ml.constrained.results, 
                                newdata = race.south, 
                                avg = TRUE, interval = "confidence")
  
  avg.pred.nonsouth.mle <- predict(ml.constrained.results, 
                                   newdata = race.nonsouth, 
                                   avg = TRUE, interval = "confidence")
  
  # A plot of a single estimate and its confidence interval
  plot(avg.pred.south.mle, labels = "South")
  
  # A  plot of the two estimates and their confidence intervals
  # use c() to combine more than one predict object for plotting
  plot(c(avg.pred.south.mle, avg.pred.nonsouth.mle), labels = c("South", "Non-South"))
  
  # The difference option can also be used to simultaneously
  # calculate separate estimates of the two sub-groups
  # and the estimated difference. This can also be plotted.
  
  avg.pred.diff.mle <- predict(ml.constrained.results, 
                               newdata = race.south, newdata.diff = race.nonsouth,
                               se.fit = TRUE, avg = TRUE)
  
  ##warning("could not do diff mle plot")
  
  ##plot(avg.pred.diff.mle)
  
  # Social desirability bias plots
  
  # Estimate logit for direct sensitive question
  
  data(mis)
  
  mis.list <- subset(mis, list.data == 1)
  
  mis.sens <- subset(mis, sens.data == 1)
  
  # Fit EM algorithm ML model
  
  fit.list <- ictreg(y ~ age + college + male + south,
                     J = 4, data = mis.list, method = "ml")
  
  # Fit logistic regression with directly-asked sensitive question
  
  fit.sens <- glm(sensitive ~ age + college + male + south, 
                  data = mis.sens, family = binomial("logit"))
  
  # Predict difference between response to sensitive item
  # under the direct and indirect questions (the list experiment).
  # This is an estimate of the revealed social desirability bias
  # of respondents. See Blair and Imai (2010).
  
  avg.pred.social.desirability <- predict(fit.list, 
                                          direct.glm = fit.sens, se.fit = TRUE)
  
  plot(avg.pred.social.desirability)
  
  
})

test_that("summary with boundary proportions for ictreg", {
  
  skip_on_cran()
  
  ceiling.results <- ictreg(y ~ age + college + male + south, treat = "treat", 
                            J = 3, data = affirm, method = "ml", fit.start = "nls",
                            ceiling = TRUE, ceiling.fit = "bayesglm",
                            ceiling.formula = ~ age + college + male + south)
  
  # Summarize fit object and generate conditional probability 
  # of ceiling liars the population proportion of ceiling liars,
  # both with standard errors.
  # Replicates Table 7 Columns 3-4 last row in Blair and Imai (2012)
  
  summary(ceiling.results, boundary.proportions = TRUE)
  
})

