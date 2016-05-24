context("PLAC fit")

test_that("PLAC() calls the right function", {
  set.seed(235711)
  dat = sim.ltrc(n = 50, Cmax = 5)$dat
  p.lbs = plr(dat)$P
  expect_equal(object = p.lbs, expected = 0.8551607, tolerance = .0005)
  expect_output(PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z1 + Z2,
                      ltrc.data = dat, td.type = "none"),
                "Calling PLAC_TI()...")

  dat = sim.ltrc(n = 50, time.dep = TRUE,
                distr.A = "binomial", p.A = 0.8, Cmax = 5)$dat
  p.lbs = plr(dat)$P
  expect_equal(object = p.lbs, expected = 0, tolerance = .0005)
  expect_output(PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
                     ltrc.data = dat, td.type = "independent",
                     td.var = "Zv", t.jump = "zeta"),
                "Calling PLAC_TD()...")
  dat = sim.ltrc(n = 50, time.dep = TRUE, Zv.depA = TRUE, Cmax = 5)$dat
  p.lbs = plr(dat)$P
  expect_equal(object = p.lbs, expected = 0.1222954, tolerance = .0005)
  expect_output(PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
                     ltrc.data = dat, td.type = "post-trunc",
                     td.var = "Zv", t.jump = "zeta"),
                "Calling PLAC_TDR()...")
})

test_that("calling cum.haz() to get the cumulative baseline function", {
  set.seed(235711)
  dat1 = sim.ltrc(n = 100, distr.T = "lnorm")$dat
  est = PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z1 + Z2,
       ltrc.data = dat1, td.type = "none")
  H = cum.haz(est, t.eval = seq(0.2, 1, 0.2))
  expect_output(str(H), "List of 5")
})
