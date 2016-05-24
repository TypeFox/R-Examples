context("Weight Calculations")

test_that("weights calculations equal", {
  XXX <- cbind(rep(1, 10), seq(0, 1, length = 10))
  AAA <- rep(c(0, 1), 5)
  fff <- c(.1, .1)
  aaa <- c(0, .9)
  GGG <- sort(rep(1:5, 2))
  mmm <- matrix(c(NA, 0.616547184192009, NA, 0.614796738538913, NA, 
                 0.613057773666495, NA, 0.611330209096351,
                 NA, 0.609613965009724), byrow = T, ncol = 2,
               dimnames = list(unique(GGG), aaa))
  
  # Checking weight calculations
  expect_equal(wght_calc(integrand = logit_integrand, 
                         allocation = aaa[2], 
                         X = XXX[1:2, ], A = AAA[1:2], 
                         fixed.effects = fff, 
                         randomization = .5), 0.462462731646267)
  
  # Checking derivative calculations
  expect_equal(wght_deriv_calc(integrand = logit_integrand, 
                               allocation = .3, 
                               X = XXX[1:2, ], A = AAA[1:2], 
                               fixed.effects = c(1, .5), 
                               random.effects = 1,
                               randomization = .5), c(-0.1277388445291851,
                                                      -0.0280714662250347,
                                                      0.0947948595448454))
  
  # Checking weight matrix calculations
  expect_equal(wght_matrix(integrand = logit_integrand,
                           allocations = aaa,
                           X = XXX, A = AAA, G = GGG,
                           fixed.effects = fff,
                           random.effects = 5,
                           randomization = .5, 
                           integrate.allocation = TRUE), mmm)
})