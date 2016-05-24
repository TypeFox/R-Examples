library(TauStar)
context("Testing the tStar function.")

# The Bergsma and Dassios (2014) 'a' function.
a = function(z) {
  sign(round(abs(z[1] - z[2]) +
               abs(z[3] - z[4]) -
               abs(z[1] - z[3]) -
               abs(z[2] - z[4]), 10))
}

# An extremely naive implementation of tStar just to check things work
# correctly in general.
tStarSlow = function(x, y, vStat = F) {
  if(length(x) != length(y) || length(x) < 4) {
    stop("Input to tStarSlow of invalid length.")
  }
  n = length(x)
  val = 0
  for(i in 1:n) {
    for(j in 1:n) {
      for(k in 1:n) {
        for(l in 1:n) {
          inds = c(i,j,k,l)
          if(length(unique(inds)) == 4 || vStat == T) {
            val = val + a(x[inds]) * a(y[inds])
          }
        }
      }
    }
  }
  if(vStat) {
    return(val / n^4)
  } else {
    return(val / (n * (n - 1) * (n - 2) * (n - 3)))
  }
}

# A distribution that is a mixture of continuous and discrete, used to check
# the tStar algorithm works on such input.
poissonGaussMix = function(n) {
  poisOrGaus = sample(c(0,1), n, replace=T)
  return(rpois(n, 5) * poisOrGaus + rnorm(n) * (1 - poisOrGaus))
}

test_that("tStar agrees with slow implementation", {
  set.seed(283721)
  m = 10
  reps = 10

  for(i in 1:reps) {
    x <- rnorm(m)
    y <- rnorm(m)
    ts = tStar(x, y)
    tvs = tStar(x, y, T)
    expect_equal(ts, tStarSlow(x, y))
    expect_equal(tvs, tStarSlow(x, y, T))
    expect_equal(ts, tStar(x, y, F, slow = T))
    expect_equal(tvs, tStar(x, y, T, slow = T))

    x <- rpois(m, 5)
    y <- rpois(m, 5)
    ts = tStar(x, y)
    tvs = tStar(x, y, T)
    expect_equal(ts, tStarSlow(x, y))
    expect_equal(tvs, tStarSlow(x, y, T))
    expect_equal(ts, tStar(x, y, F, slow = T))
    expect_equal(tvs, tStar(x, y, T, slow = T))

    x <- rnorm(m)
    y <- rpois(m, 5)
    ts = tStar(x, y)
    tvs = tStar(x, y, T)
    expect_equal(ts, tStarSlow(x, y))
    expect_equal(tvs, tStarSlow(x, y, T))
    expect_equal(ts, tStar(x, y, F, slow = T))
    expect_equal(tvs, tStar(x, y, T, slow = T))

    x <- poissonGaussMix(m)
    y <- poissonGaussMix(m)
    ts = tStar(x, y)
    tvs = tStar(x, y, T)
    expect_equal(ts, tStarSlow(x, y))
    expect_equal(tvs, tStarSlow(x, y, T))
    expect_equal(ts, tStar(x, y, F, slow = T))
    expect_equal(tvs, tStar(x, y, T, slow = T))
  }

  x = rnorm(100)
  y = rnorm(100)
  ts = tStar(x, y)
  tvs = tStar(x, y, T)
  expect_equal(ts, tStar(x, y, slow = T))
  expect_equal(tvs, tStar(x, y, T, slow = T))
  expect_true(abs(tStar(x, y, resample = T, sampleSize = 10,
                        numResamples = 10000) - ts) < 2*10^-3)
})

test_that("tStar errors on bad input", {
  x <- list(1,2,3,4)
  y <- c(1,2,3,4)
  expect_error(tStar(x, y))

  expect_error(tStar(numeric(0), numeric(0)))
  for(i in 1:3) {
    expect_error(tStar(1:i, 1:i))
  }

  expect_error(tStar(1:10, 1:9))
  expect_error(tStar(1:9, 1:10))
  expect_error(tStar(1:10, 1:10, resample=T, slow=T))
  expect_error(tStar(1:10, 1:10, resample=T, numResamples=-1))
  expect_error(tStar(1:10, 1:10, resample=T, sampleSize=-1))
  expect_error(tStar(1:10, 1:10, vStatistic=T, resample=T))
})

