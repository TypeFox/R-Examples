library(glmmsr)
context("Model fitting")

set.seed(1)
player <- 1:10
player1 <- factor(rep(2:10, 10), levels = player)
player2 <- factor(rep(1:9, 10), levels = player)

x <- rnorm(length(player))
u0 <- rnorm(length(player))

beta0 <- 1
ability0 <- beta0*x + u0
p0 <- pnorm(ability0[player1] - ability0[player2])
y <- rbinom(length(p0), 1, p0)

formula <- y ~ 0 + Sub(ability[player1] - ability[player2])
subform <- ability[player] ~ 0 + x[player] + (1 | player)
data <- list(y = y, x = x, player1 = player1, player2 = player2)

fit <- glmm(formula, subform, data = data, family = binomial,
            method = "Laplace", verbose = 0)


test_that("doesn't matter what name used for index", {
  subform_i <- ability[i] ~ 0 + x[i] + (1 | i)
  data_i <- list(x = x, player1 = player1, player2 = player2)
  fit_i <- glmm(formula, subform_i, data = data_i, family = binomial,
                method = "Laplace", verbose = 0)
  expect_equal(fit$estim, fit_i$estim)
  expect_equal(fit$Sigma, fit_i$Sigma)
})


test_that("different forms of indexing give same result", {
  player1_num <- rep(2:10, 10)
  player2_num <- rep(1:9, 10)
  data_num <- list(x = x, player1 = player1_num, player2 = player2_num)

  fit_num <- glmm(formula, subform, data = data_num, family = binomial,
                  method = "Laplace", verbose = 0)
  expect_equal(fit$estim, fit_num$estim)
  expect_equal(fit$Sigma, fit_num$Sigma)
})

test_that("OK if don't use all rows of X", {
  player1_no_1 <- factor(rep(3:11, 10), levels = 1:12)
  player2_no_1 <- factor(rep(2:10, 10), levels = 1:12)
  player1_no_1_num <- rep(3:11, 10)
  player2_no_1_num <- rep(2:10, 10)
  # has entries for 1 and 12
  data_no_1 <- list(x = c(x[1],x, x[1]),
                    player1 = player1_no_1, player2 = player2_no_1)
  data_no_1_num <- list(x = c(x[1],x, x[1]),
                    player1 = player1_no_1_num, player2 = player2_no_1_num)

  fit_no_1 <- glmm(formula, subform, data = data_no_1,
                   family = binomial, method = "Laplace", verbose = 0)
  expect_equal(fit$estim, fit_no_1$estim)
  expect_equal(fit$Sigma, fit_no_1$Sigma)

  fit_no_1_num <- glmm(formula, subform, data = data_no_1_num,
                       family = binomial, method = "Laplace", verbose = 0)
  expect_equal(fit$estim, fit_no_1_num$estim)
  expect_equal(fit$Sigma, fit_no_1_num$Sigma)
})

test_that("includes offset", {
  set.seed(1)
  offset <- rnorm(length(y))
  fit_off <- glmm(formula, subform, data = data, family = binomial,
                  offset = offset, method = "Laplace", verbose = 0)
  expect_false(identical(fit$estim, fit_off$estim))
})

test_that("random effects at observation level work OK", {
  data_re_obs <- data
  gr <- rep(1, length(y))
  gr[1:(length(y)/2)] <- 2
  gr <- as.factor(gr)
  data_re_obs$gr <- gr
  formula_re_obs <- y ~ 0 + (1 | gr) + Sub(ability[player1] - ability[player2])
  fit_re_obs <- glmm(formula_re_obs, subform, data = data_re_obs,
                     family = binomial, method = "Laplace", verbose = 0)
})

test_that("fits a two-level model correctly", {
  mod_15_glmer <- lme4::glmer(response ~ covariate + (1 | cluster),
                        data = two_level, family = binomial, nAGQ = 15)
  estim_15_glmer <- c(mod_15_glmer@theta, mod_15_glmer@beta)
  mod_15 <- glmm(response ~ covariate + (1 | cluster),
                 data = two_level, family = binomial, method = "AGQ",
                 control = list(nAGQ = 15), verbose = 0)
  estim_15 <- mod_15$estim
  expect_true(sum(abs(estim_15_glmer - estim_15)) < 0.001)

  mod_3_SR <- glmm(response ~ covariate + (1 | cluster),
                   data = two_level, family = binomial, method = "SR",
                   control = list(nSL = 3), verbose = 0)

  estim_3_SR <- mod_3_SR$estim

  expect_true(sum(abs(estim_3_SR - estim_15)) < 0.001)

  set.seed(1)
  mod_IS_100 <- glmm(response ~ covariate + (1 | cluster),
                     data = two_level, family = binomial, method = "IS",
                     control = list(nIS = 100), verbose = 0)
  estim_IS_100 <- mod_IS_100$estim
  expect_true(sum(abs(estim_IS_100 - estim_15)) < 0.2)

})


test_that("nSL = 0 gives similar result to Laplace", {
  mod_SR_0 <- glmm(response ~ covariate + (1 | cluster) + (1 | group),
                   data = three_level, family = binomial, method = "SR",
                   control = list(nSL = 0), verbose = 0)
  mod_Laplace <- glmm(response ~ covariate + (1 | cluster) + (1 | group),
                      data = three_level, family = binomial,
                      method = "Laplace", verbose = 0)

  estim_SR_0 <- mod_SR_0$estim
  estim_Laplace <- mod_Laplace$estim

  expect_true(sum(abs(estim_SR_0 - estim_Laplace)) < 1e-3)

})

test_that("factor response handled correctly", {
  two_level_factor <- two_level
  two_level_factor$response <- factor(c("N", "Y")[two_level$response + 1],
                                      levels = c("N", "Y"))
  mod_num <- glmm(response ~ covariate + (1 | cluster),
                  data = two_level, family = binomial, method = "SR",
                  control = list(nSL = 3), verbose = 0)
  mod_fac <- glmm(response ~ covariate + (1 | cluster),
                  data = two_level_factor, family = binomial, method = "SR",
                  control = list(nSL = 3), verbose = 0)

  expect_true(sum(abs(mod_num$estim - mod_fac$estim)) < 0.001)
})

test_that("warns about unused control parameters", {
  expect_warning(
    glmm(response ~ covariate + (1 | cluster),
         data = two_level, family = binomial, method = "SR",
         control = list(nSR = 3), verbose = 0),
    "unknown names"
  )
  expect_warning(
    glmm(response ~ covariate + (1 | cluster),
         data = two_level, family = binomial, method = "Laplace",
         control = list(nAGQ = 10), verbose = 0),
    "parts of control were ignored"
  )
})

test_that("uses weights", {
  two_level_double <- list(response = rep(two_level$response, 2),
                           covariate = rep(two_level$covariate, 2),
                           cluster = rep(two_level$cluster, 2))
  mod_double_Laplace <-  glmm(response ~ covariate + (1 | cluster),
                              data = two_level_double, family = binomial,
                              method = "Laplace", verbose = 0)
  mod_w2_Laplace <- glmm(response ~ covariate + (1 | cluster),
                         data = two_level, family = binomial,
                         method = "Laplace", weights = rep(2, length(two_level$response)),
                         verbose = 0)
  expect_equal(mod_double_Laplace$estim, mod_w2_Laplace$estim)

  mod_double_SR <-  glmm(response ~ covariate + (1 | cluster),
                         data = two_level_double, family = binomial,
                          method = "SR", verbose = 0)
  mod_w2_SR <- glmm(response ~ covariate + (1 | cluster),
                    data = two_level, family = binomial,
                    method = "SR", weights = rep(2, length(two_level$response)),
                    verbose = 0)
  expect_true(sum(abs(mod_double_SR$estim - mod_w2_SR$estim)) < 0.01)

  expect_error(glmm(response ~ covariate + (1 | cluster),
                    data = two_level, family = binomial,
                    method = "SR", weights = rep(2.1, length(two_level$response)),
                    verbose = 0),
               "non-integer weights")
})

test_that("checks family", {
  expect_error(glmm(response ~ covariate + (1 | cluster),
                    data = two_level, family = poisson, method = "SR", verbose = 0),
               "Only binomial family currently implemented")
  expect_error(glmm(response ~ covariate + (1 | cluster),
                    data = two_level, family = binomial(link = "cauchit"),
                    method = "SR", verbose = 0),
               "Only logit and probit links")
})

test_that("prev_fit doesn't affect the result", {
  library(BradleyTerry2)
  result <- rep(1, nrow(flatlizards$contests))
  flatlizards_glmmsr <- c(list(result = result,
                               winner = flatlizards$contests$winner,
                               loser = flatlizards$contests$loser),
                          flatlizards$predictors)
  fit_1 <- glmm(result ~ 0 + Sub(ability[winner] - ability[loser]),
                ability[liz] ~ 0 + SVL[liz] + (1 | liz),
                data = flatlizards_glmmsr, family = binomial(link = "probit"),
                method = "Laplace", verbose = 0)
  fit_2 <- glmm(result ~ 0 + Sub(ability[winner] - ability[loser]),
                ability[liz] ~ 0 + SVL[liz] + (1 | liz),
                data = flatlizards_glmmsr, family = binomial(link = "probit"),
                method = "Laplace", verbose = 0, prev_fit = fit_1)
  fit_3 <- glmm(result ~ 0 + Sub(ability[winner] - ability[loser]),
                ability[liz] ~ 0 + SVL[liz] + (1 | liz),
                data = flatlizards_glmmsr, family = binomial(link = "probit"),
                method = "Laplace", verbose = 0, prev_fit = fit_2)
  expect_true(sum(abs(fit_1$estim - fit_2$estim)) < 0.01)
  expect_true(sum(abs(fit_1$Sigma - fit_2$Sigma)) < 0.01)
  expect_true(sum(abs(fit_1$Sigma - fit_3$Sigma)) < 0.01)
  fit_SR_2 <- glmm(result ~ 0 + Sub(ability[winner] - ability[loser]),
                ability[liz] ~ 0 + SVL[liz] + (1 | liz),
                data = flatlizards_glmmsr, family = binomial(link = "probit"),
                method = "SR", control = list(nSL = 2), verbose = 0,
                prev_fit = fit_1)
  fit_SR_2_2 <- glmm(result ~ 0 + Sub(ability[winner] - ability[loser]),
                     ability[liz] ~ 0 + SVL[liz] + (1 | liz),
                     data = flatlizards_glmmsr, family = binomial(link = "probit"),
                     method = "SR", control = list(nSL = 2), verbose = 0,
                     prev_fit = fit_SR_2)
  expect_true(sum(abs(fit_SR_2$estim - fit_SR_2_2$estim)) < 0.01)
  expect_true(sum(abs(fit_SR_2$Sigma - fit_SR_2_2$Sigma)) < 0.01)
})
