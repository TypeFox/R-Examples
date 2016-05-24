
# Tests for BESTmcmc and retro power with gamma priors
#  (both work with the same BESTmcmc output object).


context("BESTmcmc&retroPower_gammaPriors")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2 <- BESTmcmc(y1, y2, priors=list(),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups and default gamma priors gives same output",  {
  expect_that(class(Bout2), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout2$mu1, 5),
      is_equivalent_to(c(4.19866, 4.10257, 4.02362, 4.61437, 4.44746, 4.35208, 5.19768, 5.24549, 4.85071)))
    expect_that(round(Bout2$mu2, 5),
      is_equivalent_to(c(2.87582, 3.26095, 3.27019, 3.38234, 3.31267, 3.25178, 3.44527, 3.42282, 3.29180)))
    expect_that(round(Bout2$nu, 5),
      is_equivalent_to(c(4.87683,  6.37871,  4.40384, 21.94421, 13.35078,  7.70271, 31.38644,  5.16176, 26.68149)))
    expect_that(round(Bout2$sigma1, 5),
      is_equivalent_to(c(0.57318, 0.90556, 0.96520, 1.03873, 1.01859, 1.12183, 0.81686, 0.85206, 0.86874)))
    expect_that(round(Bout2$sigma2, 5),
      is_equivalent_to(c(0.74529, 0.80340, 0.87125, 1.02920, 0.97957, 1.07980, 0.53805, 0.54332, 0.63603)))
    expect_that(round(attr(Bout2, "Rhat"), 5),
      is_equivalent_to(c(4.88533,  1.58644,  1.52104,  1.73854,  6.08784)))
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
      is_equivalent_to(c(0.27273, 0.09091, 0.09091, 0.18182, 0.09091,
        0.09091, 0.36364, 0.09091, 0.09091, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow2[, 2], 5),
      is_equivalent_to(c(0.04645, 0, 0, 0.00729, 0, 0, 0.10678, 0,
        0, 0, 0, 0)))
    expect_that(round(pow2[, 3], 5),
      is_equivalent_to(c(0.52243, 0.25887, 0.25887, 0.39781, 0.25887,
        0.25887, 0.63320, 0.25887, 0.25887, 0.25887, 0.25887, 0.25887)))
  }
})

Bout2a <- BESTmcmc(y1, y2,
  priors=list(muM=7:8, muSD=10:11, sigmaMode=4, sigmaSD=8, nuMean=5, nuSD=10),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups and informative gamma priors gives same output",  {
  expect_that(class(Bout2a), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2a),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
if(packageVersion("rjags") >= "4.0.0")  {
   expect_that(round(Bout2a$mu1, 5),
      is_equivalent_to(c(4.87431, 4.89931, 4.99393, 5.36281, 4.89657, 4.64511, 4.10398, 4.66114, 4.89637)))
    expect_that(round(Bout2a$mu2, 5),
      is_equivalent_to(c(3.83456, 3.44608, 3.45148, 3.06715, 3.19178, 2.87872, 3.41858, 3.48904, 3.53020)))
    expect_that(round(Bout2a$nu, 5),
      is_equivalent_to(c(9.36788, 12.99581, 11.80332, 27.46340, 14.31433, 14.18142,  5.07518,  5.34006, 9.40765)))
    expect_that(round(Bout2a$sigma1, 5),
      is_equivalent_to(c(0.68097, 0.84179, 0.75525, 1.46387, 1.05160, 0.80304, 1.49956, 1.51932, 1.55296)))
    expect_that(round(Bout2a$sigma2, 5),
      is_equivalent_to(c(0.90436, 1.29768, 0.42162, 0.34158, 0.62746, 0.48735, 0.62210, 0.53181, 0.56629)))
    expect_that(round(attr(Bout2a, "Rhat"), 5),
      is_equivalent_to(c(1.33428,  2.72412,  2.07450,  2.93150,  1.44354)))
    expect_that(attr(Bout2a, "n.eff"),
      is_equivalent_to(c(9,  9,  9,  9,  9)))
  }
  PR <- attr(Bout2a, "priors")
  expect_that(names(PR), equals(c("muM", "muSD", "sigmaMode", "sigmaSD",
    "nuMean", "nuSD")))
  expect_that(PR$muM, equals(7:8))
  expect_that(PR$muSD, equals(10:11))
})


#### One group tests
#### ===============
y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1 <- BESTmcmc(y0, priors=list(),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group and default gamma priors gives same output",  {
  expect_that(class(Bout1), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1),
    equals(c("mu", "nu", "sigma")))
if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout1$mu, 5),
      is_equivalent_to(c(1.49340, 1.26033, 1.12962, 1.60341, 1.56842, 1.55486, 1.52775, 1.44587, 1.53902)))
    expect_that(round(Bout1$nu, 5),
      is_equivalent_to(c(89.67753, 58.25130, 31.26162,  9.33254, 21.96076, 17.65331, 12.09874,  9.62349, 12.08013)))
    expect_that(round(Bout1$sigma, 5),
      is_equivalent_to(c(0.82672, 0.75241, 0.55118, 0.37248, 0.49932, 0.51823, 0.43307, 0.57858, 0.50673)))
    expect_that(round(attr(Bout1, "Rhat"), 5),
      is_equivalent_to(c(2.13151, 2.45076, 2.06028)))
    expect_that(attr(Bout1, "n.eff"),
      is_equivalent_to(c(9,  9,  9)))
  }
})

test_that("BESTpower retro with 1 group and default gamma priors gives same output",  {
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
      is_equivalent_to(c(0.72727, 0.09091, 0.09091, 0.81818, 0.09091,
        0.09091, 0.45455, 0.81818, 0.36364, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow1[, 2], 5),
      is_equivalent_to(c(0.47757, 0, 0, 0.60219, 0, 0, 0.18155, 0.60219,
        0.10678, 0, 0, 0)))
    expect_that(round(pow1[, 3], 5),
      is_equivalent_to(c(0.95355, 0.25887, 0.25887, 0.99271, 0.25887,
        0.25887, 0.73161, 0.99271, 0.63320, 0.25887, 0.25887, 0.25887)))
  }
})


y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1a <- BESTmcmc(y0,
  priors=list(muM=2, muSD=10, sigmaMode=3, sigmaSD=12, nuMean=5, nuSD=10),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group and informative gamma priors gives same output",  {
  expect_that(class(Bout1a), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1a),
    equals(c("mu", "nu", "sigma")))
if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout1a$mu, 5),
      is_equivalent_to(c(1.54808, 0.99476, 0.53802, 1.39121, 1.67442, 1.56586, 1.54330, 1.39251, 1.39831)))
    expect_that(round(Bout1a$nu, 5),
      is_equivalent_to(c(39.96915, 48.01283, 13.00077,  3.26089,  3.29042,  3.40259, 16.58101,  5.31490, 9.56252)))
    expect_that(round(Bout1a$sigma, 5),
      is_equivalent_to(c(1.98161, 1.66024, 1.59789, 0.55158, 0.29476, 0.30948, 0.55700, 0.42329, 0.39961)))
    expect_that(round(attr(Bout1a, "Rhat"), 5),
      is_equivalent_to(c(1.59382, 2.29221, 7.48744)))
    expect_that(attr(Bout1a, "n.eff"),
      is_equivalent_to(c(9,  9,  9)))
  }
  expect_that(names(attr(Bout1a, "priors")), equals(c("muM", "muSD", "sigmaMode",
    "sigmaSD", "nuMean", "nuSD")))
})
