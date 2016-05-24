data(grave)
out = spscan.test(grave, nsim = 999, alpha = .8)

test_that("check accuracy for spscan.test with SatScan for grave data", {
  expect_that(out$clusters[[1]]$locids, 
              equals(c(126, 124, 125, 132, 120, 113, 140, 141, 112, 41, 136, 30, 133, 108,
                       31, 143, 110, 139, 75)))
  expect_that(out$clusters[[1]]$coords, equals(matrix(c(10324,4389), nrow = 1)))
  expect_that(out$clusters[[1]]$pop, equals(19))
  expect_that(out$clusters[[1]]$cases, equals(11))
  expect_that(round(out$clusters[[1]]$exp, 2), equals(3.99))
  expect_that(round(out$clusters[[1]]$smr, 2), equals(2.76))
  expect_that(round(out$clusters[[1]]$rr, 2), equals(3.78))
  expect_that(round(out$clusters[[1]]$loglik, 2), equals(7.42))
  
  expect_that(out$clusters[[2]]$locids, 
              equals(c(66, 89, 71)))
  expect_that(out$clusters[[2]]$coords, equals(matrix(c(6934,6918), nrow = 1)))
  expect_that(out$clusters[[2]]$pop, equals(3))
  expect_that(out$clusters[[2]]$cases, equals(3))
  expect_that(round(out$clusters[[2]]$exp, 2), equals(0.63))
  expect_that(round(out$clusters[[2]]$smr, 2), equals(4.77))
  expect_that(round(out$clusters[[2]]$rr, 2), equals(5.19))
  expect_that(round(out$clusters[[2]]$loglik, 2), equals(4.81))
})
