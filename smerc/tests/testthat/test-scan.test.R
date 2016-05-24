set.seed(1)
coords = runif(8)

test_that("sanity checks for scan.test arguments", {
  expect_that(scan.test(coords), throws_error())
  coords = data.frame(x = runif(4), y = runif(4))
  cases = 1:3
  expect_that(scan.test(coords, cases = cases), throws_error())
  cases = 1:4
  expect_that(scan.test(coords, cases = as.factor(cases)), throws_error())
  pop = list(1:3)
  expect_that(scan.test(coords, cases = cases, pop = pop), throws_error())
  pop = list(1:4)
  expect_that(scan.test(coords, cases = cases, pop = pop), throws_error())
  pop = 1:4
  expect_that(scan.test(coords, cases = cases, pop = factor(pop)), throws_error())
  ex = 1:3
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex), throws_error())
  ex = 1:4
  alpha = -1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 1.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = c(0.1, 0.3)
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 0.1
  nsim = 0
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim), 
              throws_error())
  nsim = c(10, 20)
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim), 
              throws_error())
  nsim = 10
  ubpop = -0.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop), throws_error())
  ubpop = 1.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim,
                          ubpop = ubpop), throws_error())
  ubpop = 1.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop), throws_error())
  ubpop = c(0.1, 0.2)
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop), throws_error())
  ubpop = 0.5
  lonlat = 1:2
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, lonlat = lonlat), throws_error())
  lonlat = 1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, lonlat = lonlat), throws_error())
  lonlat = FALSE
  parallel = 1:2
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, lonlat = lonlat, parallel = parallel), throws_error())
  parallel = 1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, lonlat = lonlat, parallel = parallel), throws_error())
})

data(nydf)
out = scan.test(coords = cbind(nydf$longitude, nydf$latitude), 
                cases = floor(nydf$cases), pop = nydf$population, 
                lonlat = TRUE, nsim = 999, alpha = .50)

test_that("check accuracy for scan.test with SatScan for NY data", {
  expect_that(out$clusters[[1]]$locids , 
              equals(c(52, 50, 53, 38, 49, 48, 15, 39, 37, 1, 16, 44, 47, 40, 14, 2, 51, 13,
                       43, 45, 17, 55, 11, 3, 12, 46, 36, 35, 54, 10, 5)))
  expect_that(out$clusters[[2]]$locids , 
              equals(c(88, 87, 92, 86, 89, 91, 93, 85, 90)))
  expect_that(out$clusters[[3]]$locids , 
              equals(c(113, 117, 116, 112, 220, 118, 115, 123, 124, 111, 114, 125, 219, 122,
                       126, 119)))
  expect_that(out$clusters[[1]]$pop, equals(119050))
  expect_that(out$clusters[[1]]$cases, equals(106))
  expect_that(round(out$clusters[[1]]$exp, 2), equals(62.13))
  expect_that(round(out$clusters[[1]]$smr, 2), equals(1.71))
  expect_that(round(out$clusters[[1]]$rr, 2), equals(1.87))
  expect_that(round(out$clusters[[1]]$loglik, 2), equals(14.78))
  
  expect_that(out$clusters[[2]]$pop, equals(40696))
  expect_that(out$clusters[[2]]$cases, equals(42))
  expect_that(round(out$clusters[[2]]$exp, 2), equals(21.24))
  expect_that(round(out$clusters[[2]]$smr, 2), equals(1.98))
  expect_that(round(out$clusters[[2]]$rr, 2), equals(2.06))
  expect_that(round(out$clusters[[2]]$loglik, 2), equals(8.29))
  
  expect_that(out$clusters[[3]]$pop, equals(45667))
  expect_that(out$clusters[[3]]$cases, equals(44))
  expect_that(round(out$clusters[[3]]$exp, 2), equals(23.83))
  expect_that(round(out$clusters[[3]]$smr, 2), equals(1.85))
  expect_that(round(out$clusters[[3]]$rr, 2), equals(1.92))
  expect_that(round(out$clusters[[3]]$loglik, 2), equals(7.20))
})
