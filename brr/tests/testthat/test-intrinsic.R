context("Intrinsic loss")

test_that("Figure 10 JSPI", {
  loss <- intrinsic_phi0(phi0=0.5, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)
  expect_less_than(loss, log(1.713))
  expect_less_than(log(1.7129), loss)
  phi0 <- c(0.12, 0.1202, .994, .995)
  loss <- intrinsic_phi0(phi0=phi0, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0) 
  expect_true(loss[3] < log(10) & loss[4] > log(10))
  expect_true(loss[1] > log(10) & loss[2] < log(10))
  expect_equal(intrinsic_estimate(x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0, subdivisions=100), .4090702, tolerance=1e-4, check.attributes=FALSE)
})

test_that("Figure 15 Master", {
  x=10; y=20; S=1; T=1; a=0.5; b=0; c=0.5; d=0
  bounds <- intrinsic_bounds(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d, conf=90/100)
  loss1 <- intrinsic_phi0(bounds[1], x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d, subdivisions=100)
  loss2 <- intrinsic_phi0(bounds[2], x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d, subdivisions=100)
  expect_equal(loss1, loss2)
  expect_equal(loss1, 1.680312, tolerance=1e-5)
  estimate <- intrinsic_estimate(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d, subdivisions=100)
  expect_equal(estimate, 0.505509, tolerance=1e-4, check.attributes=FALSE)
  expect_equal(attr(estimate, "loss"), 0.4693625, tolerance=1e-5)
  mode <- spost_phi(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)$mode
  expect_less_than(mode, estimate)
})

test_that("Comparisons with simulations", {
  set.seed(666)
  phi0 <- seq(0,5,len=10)
  loss1 <- intrinsic_phi0(phi0=phi0, x=6, y=15, S=1, T=1)
  loss2 <- intrinsic_phi0_sims(phi0=phi0, x=6, y=15, S=1, T=1)
  expect_true(max(abs(loss1-loss2))<5e-2)
  #
  x <- 0; y <- 400
  loss1 <- intrinsic_phi0(phi0=phi0, x=x, y=y, S=1, T=1)
  loss2 <- intrinsic_phi0_sims(phi0=phi0, x=x, y=y, S=1, T=1)
  expect_true(max(abs(loss1-loss2))<5e-2)  
})