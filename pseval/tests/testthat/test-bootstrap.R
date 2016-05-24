library(pseval)
library(survival)

test_that("Testing bootstrap and VE estimation", {

  set.seed(52001)
  fakedata <- generate_example_data(n = 200)
  binary.ps <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP) +
    integrate_parametric(S.1 ~ BIP) + risk_binary(D = 10) + ps_estimate()

  binary.boot <- binary.ps + ps_bootstrap(n.boots = 50, progress.bar = FALSE, start = binary.ps$estimates$par)

  expect_is(binary.boot, "psdesign")
  expect_true("convergence" %in% colnames(binary.boot$bootstraps))
  expect_true(nrow(binary.boot$bootstraps) == 50)

  pw.ve <- calc_risk(binary.boot, n.samps = 100, CI.type = "pointwise")
  band.ve <- calc_risk(binary.boot, n.samps = 100, CI.type = "band")

  expect_true("Y.boot.se" %in% colnames(band.ve))
  expect_true("Y.boot.se" %in% colnames(pw.ve))

  smary <- summary(binary.boot)

  expect_equal(names(smary), c("print", "VE.estimates"))
  expect_equal(names(smary$print), c("wem.test", "boot.table"))

  expect_less_than(smary$print$wem.test$p.value, 1)
  expect_more_than(smary$print$wem.test$p.value, 0)

  expect_true(smary$print$boot.table$conv[1] == 50)
  expect_true(length(smary$VE.estimates) == 3)

})