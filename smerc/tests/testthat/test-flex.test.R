set.seed(1)
coords = runif(8)

test_that("sanity checks for flex.test arguments", {
  expect_that(flex.test(coords), throws_error())
  coords = data.frame(x = runif(4), y = runif(4))
  cases = 1:3
  expect_that(flex.test(coords, cases = cases), throws_error())
  cases = 1:4
  expect_that(flex.test(coords, cases = as.factor(cases)), throws_error())
  pop = list(1:3)
  expect_that(flex.test(coords, cases = cases, pop = pop), throws_error())
  pop = list(1:4)
  expect_that(flex.test(coords, cases = cases, pop = pop), throws_error())
  pop = 1:4
  expect_that(flex.test(coords, cases = cases, pop = factor(pop)), throws_error())
  ex = 1:3
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex), throws_error())
  ex = 1:4
  alpha = -1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 1.1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = c(0.1, 0.3)
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 0.1
  nsim = 0
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim), 
              throws_error())
  nsim = c(10, 20)
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim), 
              throws_error())
  nsim = 10
  k = 0.5
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k), throws_error())
  k = c(1, 2)
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim,
                          k = k), throws_error())
  k = 2
  lonlat = 1:2
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, lonlat = lonlat), throws_error())
  lonlat = 1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, lonlat = lonlat), throws_error())
  lonlat = FALSE
  parallel = 1:2
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, lonlat = lonlat, parallel = parallel), throws_error())
  parallel = 1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, lonlat = lonlat, parallel = parallel), throws_error())
  parallel = TRUE
  w = 1:4
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = matrix(1:4)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = diag(3)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = matrix(factor(diag(4)), nrow = 4)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = diag(4)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w, k = 10), throws_error())
})

data(nydf)
data(nyw)
out = flex.test(coords = cbind(nydf$longitude, nydf$latitude), 
                cases = floor(nydf$cases), pop = nydf$population, 
                w = nyw, k = 5,
                lonlat = TRUE, nsim = 99, alpha = .90, nreport = 5)

test_that("check accuracy for scan.test with FlexScan original for NY data", {
  expect_that(out$clusters[[1]]$locids , 
              equals(c(86, 88, 89, 92)))
  expect_that(out$clusters[[1]]$cases, equals(24))
  expect_that(round(out$clusters[[1]]$exp, 2), equals(9.24))
  expect_that(round(out$clusters[[1]]$smr, 2), equals(2.60))
  expect_that(round(out$clusters[[1]]$loglik, 2), equals(8.35))
  
  expect_that(out$clusters[[2]]$locids , 
              equals(c(1, 2, 13, 15, 49)))
  expect_that(out$clusters[[2]]$cases, equals(23))
  expect_that(round(out$clusters[[2]]$exp, 2), equals(10.03))
  expect_that(round(out$clusters[[2]]$smr, 2), equals(2.29))
  expect_that(round(out$clusters[[2]]$loglik, 2), equals(6.27))
  
  expect_that(out$clusters[[3]]$locids , 
              equals(c(37, 38, 40, 43)))
  expect_that(out$clusters[[3]]$cases, equals(21))
  expect_that(round(out$clusters[[3]]$exp, 2), equals(8.91))
  expect_that(round(out$clusters[[3]]$smr, 2), equals(2.36))
  expect_that(round(out$clusters[[3]]$loglik, 2), equals(6.05))
})
