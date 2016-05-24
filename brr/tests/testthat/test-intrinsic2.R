context("Intrinsic2 loss")

test_that("Figure 19 Master", {
  a <- 0.5; b <- 250; S <- T <- 5000; phi0 <- 0.5
  d1 <- intrinsic2_discrepancy(0.2, phi0, a=a, b=b, S=S, T=T)
  d2 <- intrinsic2_discrepancy(1, phi0, a=a, b=b, S=S, T=T)
  expect_equal(d1, 0.8563108, tol=1e-6)
  expect_equal(d2, 0.8678295, tol=1e-6)
})

test_that("Figure 20 Master", {
  a <- 4; b <- 1000; S <- T <- 5000
  x <- 4; y <- 20
  #seq(0.001, 0.9, by=0.05) %>% {plot(., intrinsic2_phi0(. , x, y, S, T, a, b), type="l")}
  l1 <- intrinsic2_phi0(0.2 , x, y, S, T, a, b, subdivisions=100)
  l2 <- intrinsic2_phi0(0.8 , x, y, S, T, a, b, subdivisions=100)
  expect_equal(l1, 0.4374621, tol=1e-5)
  expect_equal(l2, 4.276004, tol=1e-5)
  estimate <- intrinsic2_estimate(x, y, S, T, a, b, subdivisions=100)
  expect_equal(estimate, 0.2173689, tolerance=1e-4, check.attributes=FALSE)
  expect_equal(attr(estimate, "loss"), 0.4276275, tolerance=1e-5)
#   a<-2; b<-10; c<-1/2; d<-0; S<-100; T<-S; x<-0; y<-20
#   seq(0.0000000001, 0.1, by=0.001) %>% {plot(., intrinsic2_phi0(. , x, y, S, T, a, b), type="l")}
  
})

test_that("Comparisons with simulations", {
  a <- 3; b <- 14
  set.seed(666)
  phi0 <- seq(0,5,len=10)
  x <- 6; y <- 15
  loss1 <- intrinsic2_phi0(phi0=phi0, x=x, y=y, S=1, T=1, a=a, b=b)
  loss2 <- intrinsic2_phi0_sims(phi0=phi0, x=x, y=y, S=1, T=1, a=a, b=b)
  expect_true(max(abs(loss1-loss2))<1e-2)
  #
  x <- 0; y <- 150
  loss1 <- intrinsic2_phi0(phi0=phi0, x=x, y=y, S=1, T=1, a=a, b=b)
  loss2 <- intrinsic2_phi0_sims(phi0=phi0, x=x, y=y, S=1, T=1, a=a, b=b)
  expect_true(max(abs(loss1-loss2))<1e-2)  
})

# test_that("Figure 19 Master", {
#   loss <- intrinsic_phi0(phi0=0.5, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)
#   expect_less_than(loss, log(1.713))
#   expect_less_than(log(1.7129), loss)
#   phi0 <- c(0.12, 0.1202, .994, .995)
#   loss <- sapply(phi0, function(phi0) intrinsic_phi0(phi0=phi0, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)) 
#   expect_true(loss[3] < log(10) & loss[4] > log(10))
#   expect_true(loss[1] > log(10) & loss[2] < log(10))
#   expect_equal(intrinsic_estimate(x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0), .4090697, tolerance=1e-7, check.attributes=FALSE)
# })

# test_that("Figure 15 Master", {
#   x=10; y=20; S=1; T=1; a=0.5; b=0; c=0.5; d=0
#   bounds <- intrinsic_bounds(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d, conf=90/100)
#   loss1 <- intrinsic_phi0(bounds[1], x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)
#   loss2 <- intrinsic_phi0(bounds[2], x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)
#   expect_equal(loss1, loss2)
#   expect_equal(loss1, 1.680312, tolerance=1e-6)
#   estimate <- intrinsic_estimate(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)
#   expect_equal(estimate, 0.5055056, tolerance=1e-6, check.attributes=FALSE)
#   expect_equal(attr(estimate, "loss"), 0.4693617, tolerance=1e-6)
#   mode <- spost_phi(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)$mode
#   expect_less_than(mode, estimate)
# })

