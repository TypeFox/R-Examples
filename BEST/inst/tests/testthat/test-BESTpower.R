
# Tests for BESTpower



context("BESTpower_Pro")

test_that("BESTpower with 2 groups gives same output",  {
  proData <- makeData(mu1=108, sd1=17, mu2=100, sd2=15, nPerGrp=20, 
                         pcntOut=15, sdOutMult=2.0, rnd.seed=1,
                         showPlot=FALSE)
  expect_that(class(proData), equals("list"))
  expect_that(round(proData$y1, 5), 
    is_equivalent_to(c(94.91088, 109.80678, 91.06461, 135.7637, 112.48891, 91.34337, 115.39274, 120.00616, 117.01734, 100.81456, 134.22832, 113.59835, 95.00674, 65.70652, 127.115, 105.60376, 106.13228, 144.94097, 116.18516, 62.87387)))
  expect_that(round(proData$y2, 5), 
    is_equivalent_to(c(118.12939, 115.72711, 103.30548, 67.0728, 112.8777, 101.01111, 99.26143, 76.17697, 93.6024, 109.33356, 125.8485, 100.192, 108.80216, 101.05191, 77.82177, 94.71111, 95.07459, 59.30063, 130.72614, 109.97323)))
  proMCMC <- BESTmcmc(proData$y1, proData$y2, numSavedSteps=9,
      burnInSteps = 1, verbose=FALSE, rnd.seed=2)  
  if(packageVersion("rjags") >= "4.0.0")
    expect_that(round(colMeans(proMCMC), 5), 
      is_equivalent_to(c(106.39054,  99.44194,  59.25097,  21.89881,  17.49193)))
  pow2 <- BESTpower(proMCMC, N1=10, N2=10,
               ROPEm=c(-2,2) , ROPEsd=c(-2,2) , ROPEeff=c(-0.5,0.5) , 
               maxHDIWm=25.0 , maxHDIWsd=10.0 , maxHDIWeff=1.0 ,
               nRep=9, mcmcLength=1000, verbose=0, rnd.seed=3) 
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
  if(packageVersion("rjags") >= "4.0.0") {
    expect_that(round(pow2[, 1], 5), 
      is_equivalent_to(c(0.09091, 0.09091, 0.09091, 0.09091, 0.09091,
        0.09091, 0.09091, 0.09091, 0.09091, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow2[, 2], 5), 
      is_equivalent_to(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
    expect_that(round(pow2[, 3], 5), 
      is_equivalent_to(c(0.25887, 0.25887, 0.25887, 0.25887, 0.25887,
        0.25887, 0.25887, 0.25887, 0.25887, 0.25887, 0.25887, 0.25887)))
  }
})

test_that("BESTpower with 1 group gives same output",  {
  proData <- makeData(mu1=108, sd1=17, nPerGrp=20, 
                         pcntOut=15, sdOutMult=2.0, rnd.seed=4,
                         showPlot=FALSE)
  expect_that(class(proData), equals("list"))
  expect_that(names(proData), equals(c("y1", "y2")))
  expect_that(round(proData$y1, 5), 
    is_equivalent_to(c(102.74455, 86.85525, 116.85798, 110.68087, 132.43809, 112.63332, 71.39483, 93.74775, 137.89858, 135.39402, 110.0661, 98.53734, 106.22488, 97.26376, 98.92728, 101.74571, 122.5897, 140.77532, 122.08091, 61.14377)))
  expect_null(proData$y2)
  proMCMC <- BESTmcmc(proData$y1, proData$y2, numSavedSteps=9,
      burnInSteps = 1, verbose=FALSE, rnd.seed=2)  
  if(packageVersion("rjags") >= "4.0.0")
    expect_that(round(colMeans(proMCMC), 5), 
      is_equivalent_to(c(109.53919,  27.79089,  21.23356)))
  pow1 <- BESTpower(proMCMC, N1=10, N2=10,
               ROPEm=c(-2,2) , ROPEsd=c(-2,2) , ROPEeff=c(-0.5,0.5) , 
               maxHDIWm=25.0 , maxHDIWsd=10.0 , maxHDIWeff=1.0 ,
               nRep=9, mcmcLength=1000, verbose=0, rnd.seed=3) 
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
  expect_that(round(pow1[, 1], 5), 
    is_equivalent_to(c(0.90909, 0.09091, 0.09091, 0.45455, 0.90909, 0.09091, 0.09091, 0.09091, 0.90909, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow1[, 2], 5), 
    is_equivalent_to(c(0.74113, 0, 0, 0.18155, 0.74113, 0, 0, 0, 0.74113, 0, 0, 0)))
  expect_that(round(pow1[, 3], 5), 
    is_equivalent_to(c(1, 0.25887, 0.25887, 0.73161, 1, 0.25887, 0.25887, 0.25887, 1, 0.25887, 0.25887, 0.25887)))
})

