context("Score Calculations")

test_that("score calculations equal", {
  XXX <- cbind(rep(1, 10), seq(0, 1, length = 10))
  AAA <- rep(c(0, 1), 5)
  fff <- c(.1, .1)
  aaa <- c(0, .9)
  GGG <- sort(rep(1:5, 2))
  
  sss <- matrix(c(0.128180922270839, 0.0214436638223919, -0.0255666885718509,
                  0.127702581752999, 0.0497534281625182, -0.0259949349954316,
                  0.127225058823801, 0.0778510300867039, -0.0264194521475676,
                  0.126748355218138, 0.1057370155008366, -0.0268402501427183,
                  0.126272472861078, 0.1334119329472762, -0.0272573391208153),
                byrow = T, ncol = 3, dimnames = list(unique(GGG)))

  # Checking log likelihood calculations
  expect_equal(log_likelihood(integrand = logit_integrand, 
                         allocation = aaa[2], 
                         X = XXX[1:2, ], A = AAA[1:2], 
                         fixed.effects = fff, 
                         x = fff[1], pos = 1,
                         randomization = .5), -1.63675630315117)
  
  # Checking score calculations
  expect_equal(score_calc(integrand = logit_integrand, 
                               allocation = .3, 
                               X = XXX[1:2, ], A = AAA[1:2], 
                               fixed.effects = c(1, .5), 
                               random.effects = NULL,
                               randomization = .5), c(0.1032180156133930,
                                                      0.0286844121526597))
  
  # Checking score matrix calculations
  expect_equal(score_matrix(integrand = logit_integrand,
                                 allocations = aaa,
                                 X = XXX, A = AAA, G = GGG,
                                 fixed.effects = fff,
                                 random.effects = 5,
                                 randomization = .5), sss)
})


