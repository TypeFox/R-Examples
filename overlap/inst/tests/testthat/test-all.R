
# test_that code for the overlap package


# library(testthat)
# library(overlap)
# test_file("./overlap/inst/tests/test-all.R")

require(overlap) # otherwise can't find simulatedData

context("Built-in data sets")
test_that("built-in data sets are unchanged",  {
  data(simulatedData)
  expect_that(round(mean(tigerTrue), 6), equals(0.157957))
   expect_that(round(mean(pigTrue), 6), equals(0.157913))
  expect_that(round(mean(tigerObs), 6), equals(3.248677))
  expect_that(round(mean(pigObs), 6), equals(3.328342))

  data(kerinci)
  expect_that(dim(kerinci), equals(c(1098, 3)))
  expect_that(names(kerinci), equals(c("Zone", "Sps", "Time")))
  expect_that(sum(kerinci$Time), equals(540.68))
  expect_that(sum(kerinci$Zone), equals(2950))
  expect_that(nlevels(kerinci$Sps), equals(8))
  expect_that(summary(kerinci$Sps),
    is_equivalent_to(c(28, 86, 104, 273, 200, 25, 181, 201)))
  expect_that(levels(kerinci$Sps),
    equals(c("boar", "clouded", "golden", "macaque", "muntjac",
             "sambar", "tapir", "tiger")))
} )
  
context("Main computation functions")
test_that("overlapTrue gives correct answer", {
  data(simulatedData)
  expect_that(overlapTrue(tigerTrue, pigTrue), equals(0.2910917, tolerance = 1e-6))
  expect_that(overlapTrue(cbind(tigerTrue, pigTrue)), equals(0.2910917, tolerance = 1e-6))
})

test_that("densityFit gives correct answer", {
  data(simulatedData)
  expect_that(densityFit(tigerObs, c(0, pi/2, pi, 3*pi/2, 2*pi), 30), 
    equals(c(0.02440435, 0.44522913, 0.02179983, 0.50513539, 0.02440435), tolerance = 1e-7))
  expect_that(densityFit(pigObs, c(0, pi/2, pi, 3*pi/2, 2*pi), 10), 
    equals(c(7.877244e-06, 4.522317e-02, 4.622752e-01, 1.216268e-01, 7.877244e-06),
      tolerance = 1e-7))
})

test_that("getBandWidth gives correct answer", {
  data(simulatedData)
  expect_that(getBandWidth(tigerObs), equals(29.90645, tolerance = 1e-7))
  expect_that(getBandWidth(pigObs), equals(10.42065, tolerance = 1e-7))
})

test_that("overlapEst gives correct answer", {
  data(simulatedData)
  expect_that(round(overlapEst(tigerObs, pigObs), 6),
    is_equivalent_to(c(0.290862, 0.269201, 0.227500)))
  expect_that(
    round(overlapEst(tigerObs, pigObs, adjust=c(1.2, 1.5, 1)), 6), 
    is_equivalent_to(c(0.315068, 0.288488, 0.237500)))
  expect_that(
    round(overlapEst(tigerObs, pigObs, adjust=c(NA, 1, NA)), 6), 
    is_equivalent_to(c(NA_real_, 0.269201, NA_real_)))
})

context("Bootstrap functions")
test_that("resample smooth=TRUE gives correct answer", {
  data(simulatedData)
  set.seed(123)
  tigSim <- resample(tigerObs, 5, TRUE)
  expect_that(round(colMeans(tigSim), 6),
    equals(c(3.231229, 3.279921, 3.711014, 3.379395, 3.080605)))
  pigSim <- resample(pigObs, 5, TRUE)
  expect_that(round(colMeans(pigSim), 6),
    equals(c(3.267318, 3.305281, 3.434542, 3.303898, 3.296674)))
  boots <- bootEst(tigSim, pigSim)
  expect_that(round(colMeans(boots), 6),
    is_equivalent_to(c(0.365660, 0.348909, 0.328000)))
})

test_that("resample smooth=FALSE gives correct answer", {
  data(simulatedData)
  set.seed(123)
  tigSim <- resample(tigerObs, 5, FALSE)
  expect_that(round(colMeans(tigSim), 6),
    equals(c(3.277974, 3.033689, 3.332466, 3.141917, 3.645971)))
  pigSim <- resample(pigObs, 5, FALSE)
  expect_that(round(colMeans(pigSim), 6),
    equals(c(3.361474, 3.397498, 3.276355, 3.395407, 3.331279)))
  boots <- bootEst(tigSim, pigSim)
  expect_that(round(colMeans(boots), 6),
    is_equivalent_to(c(0.262617, 0.243341, 0.240000)))
})

context("Confidence intervals")
test_that("bootCI gives same results as boot.ci for common CIs", {
  require(boot)
  set.seed(123)
  dat <- runif(20)
  mean.b <- function(d,p,...) mean(d[p])
  bootout <- boot(dat, mean.b, 999)
  t0 <- bootout$t0
  bt <- as.vector(bootout$t)
  expect_that(t0, equals(mean(dat)))

  expect_that(bootCI(t0, bt)['norm',], 
    is_equivalent_to(boot.ci(bootout, 0.95, "norm")$norm[2:3]))
  expect_that(bootCI(t0, bt)['basic',], 
    is_equivalent_to(boot.ci(bootout, 0.95, "basic")$basic[4:5]))
  expect_that(bootCI(t0, bt)['perc',], 
    is_equivalent_to(boot.ci(bootout, 0.95, "perc")$perc[4:5]))

  expect_that(bootCI(t0, bt, 0.8)['norm',], 
    is_equivalent_to(boot.ci(bootout, 0.8, "norm")$norm[2:3]))
  expect_that(bootCI(t0, bt, 0.8)['basic',], 
    is_equivalent_to(boot.ci(bootout, 0.8, "basic")$basic[4:5]))
  expect_that(bootCI(t0, bt, 0.8)['perc',], 
    is_equivalent_to(boot.ci(bootout, 0.8, "perc")$perc[4:5]))

  expect_that(bootCIlogit(t0, bt)['norm',], 
    is_equivalent_to(boot.ci(bootout, 0.95, "norm", h=qlogis, hinv=plogis)$norm[2:3]))
  expect_that(bootCIlogit(t0, bt)['basic',], 
    is_equivalent_to(boot.ci(bootout, 0.95, "basic", h=qlogis, hinv=plogis)$basic[4:5]))
} )

test_that("bootCI gives correct results", {
  set.seed(123)
  dat <- runif(20)
  t0 <- sd(dat)
  bootmat <- matrix(sample(dat, 20*999, replace=TRUE), 20, 999)
  bt <- apply(bootmat, 2, sd)

  expect_that(round(bootCI(t0, bt)['norm',], 6), 
    is_equivalent_to(c(0.260553, 0.386618)))
  expect_that(round(bootCI(t0, bt)['perc',], 6), 
    is_equivalent_to(c(0.235347, 0.365483)))
  expect_that(round(bootCI(t0, bt)['basic',], 6), 
    is_equivalent_to(c(0.261460, 0.391595)))
  expect_that(round(bootCI(t0, bt)['norm0',], 6), 
    is_equivalent_to(c(0.250439, 0.376504)))
  expect_that(round(bootCI(t0, bt)['basic0',], 6), 
    is_equivalent_to(c(0.245461, 0.375597)))
} )

test_that("bootCIlogit gives correct results", {
  set.seed(123)
  dat <- runif(20)
  t0 <- sd(dat)
  bootmat <- matrix(sample(dat, 20*999, replace=TRUE), 20, 999)
  bt <- apply(bootmat, 2, sd)

  expect_that(round(bootCIlogit(t0, bt)['norm',], 6), 
    is_equivalent_to(c(0.262109, 0.394479)))
  expect_that(round(bootCIlogit(t0, bt)['perc',], 6), 
    is_equivalent_to(c(0.235347, 0.365483)))
  expect_that(round(bootCIlogit(t0, bt)['basic',], 6), 
    is_equivalent_to(c(0.265761, 0.403832)))
  expect_that(round(bootCIlogit(t0, bt)['norm0',], 6), 
    is_equivalent_to(c(0.252147, 0.382090)))
  expect_that(round(bootCIlogit(t0, bt)['basic0',], 6), 
    is_equivalent_to(c(0.244864, 0.377662)))
} )


test_that("quantileInter gives same results as norm.inter", {
  set.seed(123)
  foo <- runif(200)
  expect_that(round(quantileInter(foo), 6),
    equals(c(0.045564, 0.979369)))
  expect_that(round(quantileInter(foo[-1]), 6),
    equals(c(0.045556, 0.979822)))
  expect_that(round(quantileInter(foo, conf=0.8), 6),
    equals(c(0.141946, 0.894843)))

  diddyfoo <- runif(20)
  expect_that(quantileInter(diddyfoo),
    equals(rep(NA_real_, 2)))
} )

context("Output from plotting functions")
test_that("densityPlot gives correct output", {
  data(simulatedData)
  foo <- densityPlot(pigObs)
  expect_that(class(foo), equals("data.frame"))
  expect_that(names(foo), equals(c("x", "y"))) 
  expect_that(nrow(foo), equals(128))
  wanted <- foo$x > 0 & foo$x < 24
  expect_that(round(mean(foo$y[wanted]) * 24, 4), equals( 0.9961))

  foo <- densityPlot(tigerObs, xscale = NA, xcenter = "m", n.grid=1024)
  expect_that(class(foo), equals("data.frame"))
  expect_that(names(foo), equals(c("x", "y"))) 
  expect_that(nrow(foo), equals(1024))
  wanted <- foo$x > -pi & foo$x < pi
  expect_that(round(mean(foo$y[wanted]) * 2 * pi, 4), equals( 1.0004))
})

test_that("overlapPlot gives correct output", {
  data(simulatedData)
  foo <- overlapPlot(pigObs, tigerObs)
  expect_that(class(foo), equals("data.frame"))
  expect_that(names(foo), equals(c("x", "densityA", "densityB"))) 
  expect_that(nrow(foo), equals(128))
  wanted <- foo$x > 0 & foo$x < 24
  expect_that(round(mean(foo$densityA[wanted]) * 24, 4), equals( 1.0079))
  expect_that(round(mean(foo$densityB[wanted]) * 24, 4), equals( 1.0067))

  foo <- overlapPlot(pigObs, tigerObs, xscale = NA, xcenter = "m", n.grid=1024)
  expect_that(class(foo), equals("data.frame"))
  expect_that(names(foo), equals(c("x", "densityA", "densityB"))) 
  expect_that(nrow(foo), equals(1024))
  wanted <- foo$x > -pi & foo$x < pi
  expect_that(round(mean(foo$densityA[wanted]) * 2 * pi, 4), equals(0.9981))
  expect_that(round(mean(foo$densityB[wanted]) * 2 * pi, 4), equals(1.0008))
})

graphics.off()


