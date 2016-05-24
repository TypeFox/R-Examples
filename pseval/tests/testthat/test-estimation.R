library(pseval)
library(survival)

test_that("Testing all combinations of integration and risk models", {

  set.seed(500)
  fakedata <- generate_example_data(n = 500)

  binary.ps <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)

  expect_is(psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP, CPV = CPV) + integrate_parametric(S.1 ~ BIP) +
    risk_binary(D = 10, risk = risk.logit) +
    ps_estimate(), "psdesign")

  binfit1 <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP, CPV = CPV) + integrate_parametric(S.1 ~ BIP) +
    risk_binary(D = 10, risk = risk.logit) +
    ps_estimate() #+ ps_bootstrap()

  #stg <- calc_STG(binfit1, permute.times = 1000)

  #plot(binfit1, contrast = "RD")
  #stg
  #hist(stg$permutation$permuted.stats)

  expect_is(calc_STG(binfit1, permute = FALSE)$obsSTG, "numeric")
  expect_true(calc_STG(binfit1, permute = FALSE)$obsSTG > 0)

  expect_is(psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP, BSM = BSM) + integrate_parametric(S.1 ~ BIP) +
              risk_binary(D = 10, risk = risk.logit) +
              ps_estimate(), "psdesign")

  fakedata$weights <- runif(nrow(fakedata))

  expect_is(psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP, CPV = CPV, BSM = BSM) + integrate_parametric(S.1 ~ BIP) +
              risk_binary(D = 10, risk = risk.logit) +
              ps_estimate(), "psdesign")

  expect_is(psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP, CPV = CPV, BSM = BSM, weights = weights) + integrate_parametric(S.1 ~ BIP) +
              risk_binary(D = 10, risk = risk.logit) +
              ps_estimate(), "psdesign")

  expect_is(binary.ps, "psdesign")

  expect_is(binary.ps +
    integrate_parametric(S.1 ~ BIP) +
    risk_binary(model = Y ~ S.1 * Z, D = 10, risk = risk.logit) +
    ps_estimate(method = "BFGS"), "psdesign")

  expect_is(binary.ps +
              integrate_parametric(S.1 ~ BIP) +
              risk_binary(model = Y ~ S.1 * Z, D = 10, risk = risk.probit) +
              ps_estimate(method = "BFGS"), "psdesign")

  expect_is(binary.ps +
              integrate_bivnorm("S.1") +
              risk_binary(model = Y ~ S.1 * Z, D = 10, risk = risk.logit) +
              ps_estimate(method = "BFGS"), "psdesign")

  expect_is(binary.ps +
              integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) +
              risk_binary(model = Y ~ S.1 * Z, D = 10, risk = risk.logit) +
              ps_estimate(method = "BFGS"), "psdesign")

  expect_error(binary.ps +
              integrate_nonparametric(S.1 ~ BIP))


  cat.ps <- psdesign(fakedata, Z = Z, Y = Y.obs,
                         S = S.obs.cat, BIP = BIP.cat)

  cat.ps.num <- psdesign(fakedata, Z = Z, Y = Y.obs,
                     S = as.numeric(S.obs.cat), BIP = as.numeric(BIP.cat))

  expect_is(cat.ps.num + integrate_nonparametric(S.1 ~ BIP) + risk_binary(D = 10) + ps_estimate(method = "pseudo-score"), "psdesign")

  ## categorical W with continuous S

  catw.ps <- psdesign(fakedata, Z = Z, Y = Y.obs,
                     S = S.obs, BIP = BIP.cat)
  expect_is(catw.ps + integrate_parametric(S.1 ~ BIP) + risk_binary(D = 10) + ps_estimate(method = "pseudo-score"), "psdesign")

  expect_is(cat.ps +
    integrate_nonparametric(formula = S.1 ~ BIP) +
    risk_binary(Y ~ S.1 * Z, D = 10, risk = risk.logit) +
      ps_estimate(method = "pseudo-score"), "psdesign")

  expect_is(cat.ps +
              integrate_nonparametric(formula = S.1 ~ BIP) +
              risk_binary(Y ~ S.1 * Z, D = 10, risk = risk.logit) +
              ps_estimate(method = "BFGS"), "psdesign")

  expect_error(cat.ps +
              integrate_parametric(formula = S.1 ~ BIP))

  expect_error(cat.ps +
                 integrate_bivnorm("S.1"))

  expect_error(cat.ps +
                 integrate_semiparametric(formula.location = S.1 ~ BIP, formula.scale = S.1 ~ 1))

  expect_warning(psdesign(fakedata, Z = Z, Y = Surv(time.obs, event.obs),
                      S = S.obs, BIP = BIP))

  surv.ps <- psdesign(fakedata, Z = Z, Y = Surv(time.obs, event.obs),
                      S = S.obs, BIP = BIP, tau = 0)

  expect_is(surv.ps +
              integrate_parametric(S.1 ~ BIP) +
              risk_exponential(D = 10) + ps_estimate(), "psdesign")

  expect_is(surv.ps +
              integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) +
              risk_exponential(D = 10) + ps_estimate(), "psdesign")

  expect_is(surv.ps +
              integrate_bivnorm() +
              risk_exponential(D = 10) + ps_estimate(), "psdesign")

  expect_error(surv.ps + integrate_nonparametric(S.1 ~ BIP))

  expect_is(surv.ps +
              integrate_parametric(S.1 ~ BIP) +
              risk_weibull(D = 10) + ps_estimate(), "psdesign")

  expect_is(surv.ps +
              integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) +
              risk_weibull(D = 10) + ps_estimate(), "psdesign")

  expect_is(surv.ps +
              integrate_bivnorm() +
              risk_weibull(D = 10) + ps_estimate(), "psdesign")

  expect_error(surv.ps + integrate_parametric(S.1 ~ BIP) + risk_weibull(D = 10) + ps_estimate(method = "pseudo-score"))


  ## count data

  fakedata.count <- fakedata
  fakedata.count$Yct <- with(fakedata.count, floor(time.obs * 100))

  count.ps <- psdesign(data = fakedata.count, Z = Z, Y = Yct, S = S.obs, BIP = BIP, timeon = time.obs)

  expect_is(count.ps + integrate_parametric(S.1 ~ BIP) +
              risk_poisson(D = 10) +
              ps_estimate(), "psdesign")

  ## include offset

  expect_is(count.ps + integrate_parametric(S.1 ~ BIP) +
              risk_poisson(model = Y ~ S.1 * Z + offset(log(timeon)), D = 10) +
              ps_estimate(), "psdesign")



})