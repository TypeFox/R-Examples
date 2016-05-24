context("Excursions")


test_that("Excursions, alpha = 1, type = >", {
  data <- integration.testdata1()
  res <- excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',
                    seed = data$seed, max.threads = 1)
  r <- c(2.463175e-15, 1.030394e-09, 7.734328e-06, 0.002534894, 0.07579555,
         0.4188056, 0.8192485, 0.9746894, 0.9984611, 0.9999625, 0.9999997)
})

test_that("Excursions, alpha = 1, type = <", {
  data <- integration.testdata1()
  res = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='<',
                   seed = data$seed, max.threads = 1)
  r = c(0.9999997, 0.9999619, 0.9984783, 0.9746628, 0.819054,
        0.4196815, 0.07603522,0.002548885,7.762954e-06,1.033062e-09,2.47039e-15)
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 1, type = =", {
  data <- integration.testdata1()
  res = excursions(alpha=1, u=0, mu=data$mu+0.1, Q=data$Q, type="=",
                   seed = data$seed, max.threads = 1)
  r = c(7.381175e-07, 8.137438e-05, 0.003200957, 0.05128127, 0.3331441,
        0.6420603, 0.1815013, 0.02201733, 0.001153959, 2.520208e-05,
        1.945911e-07)

  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 1, type = !=", {
  data <- integration.testdata1()
  res = excursions(alpha=1, u=0, mu=data$mu+0.1, Q=data$Q, type='!=',
                   seed = data$seed, max.threads = 1)
  r = c(0.9999993, 0.9999186, 0.996799, 0.9487187, 0.6668559, 0.3579397,
        0.8184987, 0.9779827, 0.998846, 0.9999748, 0.9999998)

  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = >", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='>',
                    seed = data$seed, max.threads = 1)
  r = c(0,0,0,0,0,0,0,0.9801319,0.9988921,0.9999753,0.9999998)
  res$F[is.na(res$F)] = 0
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = <", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='<',
                    seed = data$seed, max.threads = 1)
  r = c(0.9999429,0.9999434,0.9978837,0.9679513,0,0,0,0,0,0,0)
  res$F[is.na(res$F)] = 0
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = =", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='=',
                    seed = data$seed, max.threads = 1)
  r = c(7.381175e-07, 8.137438e-05, 0.003200957, 0.05128127, 1, 1,
        1, 0.02201733, 0.001153959, 2.520208e-05, 1.945911e-07)
  res$F[is.na(res$F)] = 1
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, alpha = 0.1, type = !=", {
  data <- integration.testdata1()
  res = excursions(alpha=0.1, u=0, mu=data$mu+0.1, Q=data$Q, type='!=',
                   seed = data$seed, max.threads = 1)
  r = c(0.9999993, 0.9999186, 0.996799, 0.9487187, 0, 0, 0,
        0.9779827,0.998846,0.9999748,0.9999998)
  res$F[is.na(res$F)] = 0
  expect_equal(res$F,r,tolerance=1e-7)
})

test_that("Excursions, move u to mu", {
  data <- integration.testdata1()

  res = excursions(alpha=0.1, u=1, mu=data$mu, Q=data$Q, type='>',
                   seed = data$seed, max.threads = 1)
  res2 = excursions(alpha=0.1, u=0, mu=data$mu-1, Q=data$Q, type='>',
                    seed = data$seed, max.threads = 1)
  res$F[is.na(res$F)] = 0
  res2$F[is.na(res2$F)] = 0
  expect_equal(res$F,res2$F,tolerance=1e-7)
})

test_that("Excursions, input variances", {
  data <- integration.testdata1()

  vars = diag(solve(data$Q))
  res1 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',
                    seed = data$seed, vars = vars, max.threads = 1)
  res2 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',
                    seed = data$seed, max.threads = 1)
  expect_equal(res1$F,res2$F,tolerance=1e-7)
})

test_that("Excursions, ind argument order", {
  data <- integration.testdata1()

  vars = diag(solve(data$Q))

  ind1 = c(1,2,3,4)
  ind2 = c(4,3,2,1)
  ind3 = rep(FALSE,length(data$mu))
  ind3[1:4] = TRUE
  res1 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',
                    seed = data$seed, ind = ind1, max.threads = 1)
  res2 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',
                    seed = data$seed, ind = ind2, max.threads = 1)
  res3 = excursions(alpha=1, u=0, mu=data$mu, Q=data$Q, type='>',
                    seed = data$seed, ind = ind3, max.threads = 1)

  expect_equal(res1$F,res2$F,tolerance=1e-7)
  expect_equal(res2$F,res3$F,tolerance=1e-7)
})

#Tests to add:

#Test that Q.chol and Q gives the same result

#Test the max.size argument

#Test reo

#Test rho

#Test max.threads

#Test QC method

#Test ind argument


#res1 = excursions2::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=', seed = seed, max.threads = 1)
#res2 = excursions::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=', max.threads = 1)

#plot(res1$F)
#lines(res2$F, col=2)
