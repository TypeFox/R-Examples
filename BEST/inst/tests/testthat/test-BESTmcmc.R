
# Tests for BESTmcmc and retro power with priors=NULL
#  (both work with the same BESTmcmc output object).


context("BESTmcmc&retroPower")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2 <- BESTmcmc(y1, y2, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups gives same output",  {
  expect_that(class(Bout2), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout2$mu1, 5),
      is_equivalent_to(c(4.12535, 4.95591, 5.34822, 4.98114, 4.88526, 4.25582, 4.79039, 4.71482, 4.70673)))
    expect_that(round(Bout2$mu2, 5),
      is_equivalent_to(c(2.82018, 3.47593, 2.86840, 3.19374, 3.28132, 3.17236, 3.56339, 3.49283, 3.44266)))
    expect_that(round(Bout2$nu, 5),
      is_equivalent_to(c(42.28513, 50.22671, 13.62337, 13.83349,  5.33221,  5.81619, 42.31951, 61.54992, 56.15759)))
    expect_that(round(Bout2$sigma1, 5),
      is_equivalent_to(c(0.93800, 0.85199, 1.07788, 0.62254, 0.76634, 0.60601, 0.43419, 0.40904, 0.75910)))
    expect_that(round(Bout2$sigma2, 5),
      is_equivalent_to(c(0.84074, 0.54759, 0.75758, 1.35774, 2.27456, 1.85109, 0.47483, 0.73277, 0.47497)))
    expect_that(round(attr(Bout2, "Rhat"), 5),
      is_equivalent_to(c(0.95004,  1.74919,  2.71178,  2.41296,  3.62399)))
    expect_that(attr(Bout2, "n.eff"),
      is_equivalent_to(c(9,  9,  9,  9,  9)))
  }
})

test_that("BESTpower retro with 2 groups gives same output",  {
  pow2 <- BESTpower(Bout2,
    ROPEm=c(-0.1,0.1), ROPEsd=c(-2,2), ROPEeff=c(-0.5,0.5),
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456)
  expect_that(class(pow2), equals("matrix"))
  expect_that(colnames(pow2),
    equals(c("mean", "CrIlo", "CrIhi")))
  expect_that(rownames(pow2),
    equals(c("  mean:   HDI > ROPE", "  mean:   HDI < ROPE",
      "  mean:  HDI in ROPE", "  mean: HDI width ok",
      "    sd:   HDI > ROPE", "    sd:   HDI < ROPE",
      "    sd:  HDI in ROPE", "    sd: HDI width ok",
      "effect:   HDI > ROPE", "effect:   HDI < ROPE",
      "effect:  HDI in ROPE", "effect: HDI width ok")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(pow2[, 1], 5),
      is_equivalent_to(c(0.36364, 0.09091, 0.09091, 0.36364, 0.09091, 0.09091,
        0.63636, 0.27273, 0.18182, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow2[, 2], 5),
      is_equivalent_to(c(0.10678, 0, 0, 0.10678, 0, 0, 0.36680,
        0.04645, 0.00729, 0, 0, 0)))
    expect_that(round(pow2[, 3], 5),
      is_equivalent_to(c(0.63320, 0.25887, 0.25887, 0.6332, 0.25887,
        0.25887, 0.89322, 0.52243, 0.39781, 0.25887, 0.25887, 0.25887)))
  }
})

y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1 <- BESTmcmc(y0, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group gives same output",  {
  expect_that(class(Bout1), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1),
    equals(c("mu", "nu", "sigma")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout1$mu, 5),
      is_equivalent_to(c(1.92146, 1.88927, 0.96934, 1.31110, 1.42774, 1.31009, 1.49958, 1.66729, 1.66788)))
    expect_that(round(Bout1$nu, 5),
      is_equivalent_to(c(23.60838,  4.10204, 31.07991, 23.43798, 15.86208, 43.93298, 33.19201,  8.74280, 7.68395)))
    expect_that(round(Bout1$sigma, 5),
      is_equivalent_to(c(0.80861, 0.72334, 1.15983, 0.57750, 0.85844, 0.57488, 0.29597, 0.52261, 0.54746)))
    expect_that(round(attr(Bout1, "Rhat"), 5),
      is_equivalent_to(c(1.13417, 0.96993, 1.95318)))
    expect_that(attr(Bout1, "n.eff"),
      is_equivalent_to(c(9,  9,  9)))
  }
})

test_that("BESTpower retro with 1 group gives same output",  {
  pow1 <- BESTpower(Bout1,
    ROPEm=c(-0.5,0.5), ROPEsd=c(-1,1), ROPEeff=c(-1,1),
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456)
  expect_that(class(pow1), equals("matrix"))
  expect_that(colnames(pow1),
    equals(c("mean", "CrIlo", "CrIhi")))
  expect_that(rownames(pow1),
    equals(c("  mean:   HDI > ROPE", "  mean:   HDI < ROPE",
      "  mean:  HDI in ROPE", "  mean: HDI width ok",
      "    sd:   HDI > ROPE", "    sd:   HDI < ROPE",
      "    sd:  HDI in ROPE", "    sd: HDI width ok",
      "effect:   HDI > ROPE", "effect:   HDI < ROPE",
      "effect:  HDI in ROPE", "effect: HDI width ok")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(pow1[, 1], 5),
      is_equivalent_to(c(0.63636, 0.09091, 0.09091, 0.90909, 0.09091,
        0.09091, 0.36364, 0.90909, 0.36364, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow1[, 2], 5),
      is_equivalent_to(c(0.36680, 0, 0, 0.74113, 0, 0, 0.10678, 0.74113,
        0.10678, 0, 0, 0)))
    expect_that(round(pow1[, 3], 5),
      is_equivalent_to(c(0.89322, 0.25887, 0.25887, 1, 0.25887, 0.25887,
        0.63320, 1, 0.63320, 0.25887, 0.25887, 0.25887)))
  }
})
