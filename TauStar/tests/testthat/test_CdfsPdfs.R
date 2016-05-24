library(TauStar)
context("Testing the CDFs and PDFs.")

simulateDensityAndTailProbs = function(rDist, n, x = seq(-1.1, 4, length = 20)) {
  vals = rDist(n)
  lowerTailProbs = numeric(length(x))
  for (i in 1:length(x)) {
    lowerTailProbs[i] = mean(vals <= x[i])
  }
  return(list(density = density(vals, from = -1.1, to = 4),
              lowerTailProbs = lowerTailProbs,
              x = x))
}

read.matrix = function(file) {
  return(read.table(file, sep = ",", as.is = T))
}

test_that("pHoeffInd function returns correct values", {
  # Comparing against the values given by the BKR paper
  bkrToHoeff = function(x) {
    return(as.numeric(pHoeffInd(72/pi^4 * x - 1, 10^-6)))
  }

  expect_equal(.00000, bkrToHoeff(.3), tolerance = 10^-5)
  expect_true(abs(.04867 - bkrToHoeff(.6)) < 10^-5)
  expect_equal(.29652, bkrToHoeff(.9), tolerance = 10^-5)
  expect_equal(.54354, bkrToHoeff(1.2), tolerance = 10^-5)
  expect_equal(.70763, bkrToHoeff(1.5), tolerance = 10^-5)
  expect_equal(.80922, bkrToHoeff(1.8), tolerance = 10^-5)
  expect_equal(.86406, bkrToHoeff(2.05), tolerance = 10^-5)
  expect_equal(.98546, bkrToHoeff(3.9), tolerance = 10^-5)
  expect_equal(.98965, bkrToHoeff(4.2), tolerance = 10^-5)
  expect_equal(.99261, bkrToHoeff(4.5), tolerance = 10^-5)
  expect_equal(.99471, bkrToHoeff(4.8), tolerance = 10^-5)
  expect_equal(.99994, bkrToHoeff(9), tolerance = 10^-5)

  expect_equal(.95, bkrToHoeff(2.844), tolerance = 10^-4)
  expect_equal(.99, bkrToHoeff(4.230), tolerance = 10^-5)
  expect_equal(.995, bkrToHoeff(4.851), tolerance = 10^-5)

  # Just some sanity checks
  expect_true(abs(pHoeffInd(-2) - 0) < 10^-5)
  expect_true(abs(pHoeffInd(10) - 1) < 10^-6)
})

test_that("dHoeffInd function returns correct values", {
  # Comparing against the values given by the BKR paper
  empDens = read.matrix("contDensityData.csv")
  set.seed(2341)
  inds = sample(nrow(empDens), 5)
  for (i in 1:length(inds)) {
    expect_true(abs(empDens[inds[i],2] -
                      dHoeffInd(empDens[inds[i],1])) < 5*10^-3)
  }
})

test_that("pDisHoeffInd function returns correct values", {
  # Bernoulli case
  p = 2/3
  q = 1/4
  coef = 4 * p * q * (1 - p) * (1 - q)
  for (x in seq(-.2, 4, length = 10)) {
    expect_equal(as.numeric(pDisHoeffInd(x, c(p, 1 - p), c(q, 1 - q), 10^-4)),
                 pchisq(x / coef + 1, 1), tolerance = 10^-3)
  }
  lowerTailProbsMat = read.matrix("disLowerTailProbs.csv")
  p = c(0.19556498, 0.2431002, 0.10629372, 0.21154617, 0.03182509, 0.21166985)
  q = c(0.34460139, 0.29669263, 0.11664926, 0.24205672)
  for (i in 1:nrow(lowerTailProbsMat)) {
    expect_true(abs(lowerTailProbsMat[i,2] -
                 pDisHoeffInd(lowerTailProbsMat[i,1], p, q)) < 10^-4)
  }
})

test_that("dDisHoeffInd function returns correct values", {
  # Bernoulli case
  p = 2/3
  q = 1/4
  coef = 4 * p * q * (1 - p) * (1 - q)
  for (x in seq(-.2, 4, length = 10)) {
    expect_true(abs(as.numeric(dDisHoeffInd(x, c(p, 1 - p), c(q, 1 - q))) -
                 1/coef * dchisq(x / coef + 1, 1)) < 5*10^-3)
  }

  empDens = read.matrix("disDensityData.csv")
  p = c(0.19556498, 0.2431002, 0.10629372, 0.21154617, 0.03182509, 0.21166985)
  q = c(0.34460139, 0.29669263, 0.11664926, 0.24205672)
  set.seed(23784)
  inds = sample(nrow(empDens), 5)
  for (i in 1:length(inds)) {
    expect_true(abs(empDens[inds[i],2] -
                       dDisHoeffInd(empDens[inds[i],1], p, q)) < 5*10^-3)
  }
})

test_that("pMixHoeffInd function returns correct values", {
  lowerTailProbsMat = read.matrix("mixLowerTailProbs.csv")
  p = c(0.19556498, 0.2431002, 0.10629372, 0.21154617, 0.03182509, 0.21166985)
  set.seed(29923)
  for (i in 1:nrow(lowerTailProbsMat)) {
    expect_true(abs(lowerTailProbsMat[i,2] -
                      pMixHoeffInd(lowerTailProbsMat[i,1], p)) < 10^-3)
  }
})

test_that("dMixHoeffInd function returns correct values", {
  empDens = read.matrix("mixDensityData.csv")
  p = c(0.19556498, 0.2431002, 0.10629372, 0.21154617, 0.03182509, 0.21166985)
  set.seed(434)
  inds = sample(nrow(empDens), 5)
  for (i in 1:length(inds)) {
    expect_true(abs(empDens[inds[i],2] -
                      dMixHoeffInd(empDens[inds[i],1], p)) < 10^-2)
  }
})
