context("Testing {dprw}LambertW functions\n")
set.seed(10)
cauchy.samples <- rcauchy(1000)

##

beta.list <- list("normal" = c(1, 2),
                  "t" = c(1, 10, 5),
                  "exp" = 10,
                  "gamma" = c(2, 3),
                  "unif" = c(-1, 1),
                  "chisq" = 5)

dist.names <- names(beta.list)
beta.list <- lapply(names(beta.list),
                    function(nn) {
                      x <- beta.list[[nn]]
                      names(x) <- get_beta_names(nn)
                      return(x)
                    })
names(beta.list) <- dist.names

heavy.theta.list <- lapply(beta.list,
                           function(x) {
                             return(list(beta = x, delta = 0.1))
                           })
names(heavy.theta.list) <- names(beta.list)

for (nn in names(heavy.theta.list)) {
  info.txt <- paste0("Testing '", nn, "' distribution\n")
  
  context(info.txt)
  theta.tmp <- heavy.theta.list[[nn]]
  tau.tmp <- theta2tau(theta.tmp, distname = nn)
  
  theta.zero.tmp <- theta.tmp
  theta.zero.tmp[["delta"]] <- 0
  tau.zero.tmp <- theta2tau(theta.zero.tmp, distname = nn)
  
  auxD <- function(x) {
    dLambertW(x, theta = theta.tmp, distname = nn)
  }
  auxP <- function(q) {
    pLambertW(q, theta = theta.tmp, distname = nn)
  }
  auxR <- function(n) {
    rLambertW(n = n, theta = theta.tmp, distname = nn)
  }
  auxQ <- function(x) {
    qLambertW(x, theta = theta.tmp, distname = nn)
  }

  dist.family <- get_distname_family(nn)
  ib <- c(-Inf, Inf)
  if (nn %in% c("unif", "beta")) {
    ib <- theta.tmp$beta
  }
  support.dist <- get_support(tau.tmp, 
                              is.non.negative = dist.family$is.non.negative,
                              input.bounds = ib)

  test_that("dLamberW is a density", {
    area.curve <- integrate(auxD, support.dist[1], support.dist[2])$value
    expect_equal(area.curve, 1,
                 info = info.txt, tol = 1e-4)
    
    expect_true(all(auxD(cauchy.samples) >= 0))
  })
  
  
  test_that("pLamberW is a cdf", {

    expect_equal(auxP(support.dist[1]), 0,
                 info = info.txt, tol = 1e-5)
    expect_equal(auxP(support.dist[2]), 1,
                      info = info.txt, tol = 1e-5)
    
    expect_true(all(auxP(cauchy.samples) >= 0))
    expect_true(all(auxP(cauchy.samples) <= 1))
    
    # increasing function
    expect_true(all(diff(auxP(sort(cauchy.samples))) >= 0))
    
    # int_a^b pdf(x) dx = cdf(b) - cdf(a)
    expect_equal(integrate(auxD, 0, 1)$value,
                 auxP(1) - auxP(0))
  })
  
  test_that("qLambertW is a quantile function", {
    
    expect_true(auxQ(0) == support.dist[1],
                info = paste0(info.txt, ": at 0 it equals lower bound."))
    expect_true(auxQ(1) == support.dist[2],
                info = paste0(info.txt, ": at 1 it equals upper bound."))

    # qLambertW is the inverse of pLambertW
    samples.from.dist <- auxR(n = 1e2)
    expect_equal(auxQ(auxP(samples.from.dist)), samples.from.dist)
    
    # increasing function
    expect_true(all(diff(auxQ(seq(0, 1, by = 0.1))) >= 0))
    
    dist.family <- get_distname_family(nn)
    # 0 quantile must be zero for non-negative RVs
    if (dist.family$is.non.negative) {
      expect_equal(auxQ(0), 0)
    }
  })
  
  
  test_that("rLamberW is a random number generator", {

    expect_true(is.numeric(auxR(n = 100)),
                 info = info.txt)
    expect_equal(length(auxR(n = 100)), 100,
                 info = info.txt)
    
    # KDE estimated is close to true density
    samples.from.dist <- auxR(n = 1e3)
    kde.est <- density(samples.from.dist)
    expect_gt(cor(kde.est$y, auxD(kde.est$x)), 0.9)
  })
}


test_that("mLambertW only works for normal so far", {
  expect_error(mLambertW(theta = list(beta = 1), distname = "exp"))
})
test_that("mLambertW has correct values for Normal and delta = 0 or gamma = 0", {
  
  mom.lamw.h <- mLambertW(theta = list(beta = c(2, 3)),
                          distname = "normal")
  expect_equal(mom.lamw.h$mean, 2)
  expect_equal(mom.lamw.h$sd, 3)
  expect_equal(mom.lamw.h$skewness, 0)
  expect_equal(mom.lamw.h$kurtosis, 3)
  
})

test_that("mLambertW has correct values for Normal and delta > 0", {
  
  delta.v <- seq(0, 1.2, by = 0.2)
  for (dd in delta.v) {
    msg <- paste0("delta=", dd)
    mom.lamw.h <- mLambertW(theta = list(beta = c(2, 3), delta = dd),
                            distname = "normal")
    if (dd < 1) {
      expect_equal(mom.lamw.h$mean, 2, 
                   info = paste0(msg, ": mean"))
    } else {
      expect_true(is.na(mom.lamw.h$mean))
    }
    
    if (dd > 0) {
      # testthat 0.11.* removed the 'info' field for expect_gt :( 
      # will include again once fixed
      expect_gt(mom.lamw.h$sd, 3) # info = paste0(msg, ": sd"))
    } else if (dd > 0.5) {
      expect_true(is.na(mom.lamw.h$sd))
    }
    
    if (dd < 1/3) {
      expect_equal(mom.lamw.h$skewness, 0, 
                   info = paste0(msg, ": skewness"))
    } else {
      expect_true(is.na(mom.lamw.h$skewness))
    }
    
    if (dd > 0) {
      # info = paste0(msg, ": kurtosis")
      expect_gt(mom.lamw.h$kurtosis, 3)
    } else if (dd > 1/4) {
      expect_true(is.na(mom.lamw.h$mean))
    }
  }
})