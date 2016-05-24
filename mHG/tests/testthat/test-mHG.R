library(testthat)
library(mHG)

context("Utils Calculation")

test_that("d_ratio and v_ratio consistent with R", {
  N.ratios <- 100
  N <- 100
  B <- 40
  
  b <- numeric(N.ratios)
  
  d_ratio.exp <- numeric(N.ratios)
  v_ratio.exp <- numeric(N.ratios)
  d_ratio <- numeric(N.ratios)
  v_ratio <- numeric(N.ratios)
  
  n <- sample(1:N, N.ratios, replace = T)
  for (i in seq(N.ratios)) {
    if (n[i] == N) {
      b[i] <- B
    } else {
      b[i] <- sample(max(0, B + n[i] - N):min(n[i], B), 1)
    }
    HG <- dhyper(b[i], B, N - B, n[i])
    d_ratio.exp[i] <- HG / dhyper(b[i] - 1, B, N - B, n[i] - 1)
    v_ratio.exp[i] <- HG / dhyper(b[i], B, N - B, n[i] - 1)
    d_ratio[i] <- d_ratio(n[i],b[i],N,B)
    v_ratio[i] <- v_ratio(n[i],b[i],N,B)
  }
  
  expect_equal(d_ratio, d_ratio.exp, info = "d_ratio", tolerance = 1e-4)
  expect_equal(v_ratio, v_ratio.exp, info = "v_ratio", tolerance = 1e-4)
})

context("Statistic Calculation")

test_that("HG row calculation is consistent with R", {
  N.rows <- 100
  N <- 100
  B <- 40

  m <- numeric(N.rows)
  n <- numeric(N.rows)
  b_n <- numeric(N.rows)
  HG_row_m <- matrix(0, nrow = N.rows, ncol = (B + 1))
  
  HG_row_n.exp <- matrix(0, nrow = N.rows, ncol = (B + 1))
  HG_row_n.iter <- matrix(nrow = N.rows, ncol = (B + 1))
  HG_row_n.recur <- matrix(nrow = N.rows, ncol = (B + 1))
  HG_row_n.either <- matrix(nrow = N.rows, ncol = (B + 1))
  
  m <- sample(1:(N - 1), N.rows, replace = T)
  for (i in seq(N.rows)) {
    b_n[i] <- sample(max(0, B + m[i] - N):min(m[i], (B - 1)), 1) + 1
    if (m[i] == (N-1)) {
      n[i] <- N
    } else {
      n[i] <- sample((m[i]+1):N, 1)
    }
    HG_row_m[i, 1:b_n[i]] <- dhyper(0:(b_n[i] - 1), B, N - B, m[i])
    HG_row_n.exp[i, 1:(b_n[i] + 1)] <- dhyper(0:b_n[i], B, N - B, n[i])
    HG_row_n.iter[i,] <- HG_row_n.calc.iter(HG_row_m[i,], m[i], n[i], b_n[i], N, B)
    HG_row_n.recur[i,] <- HG_row_n.calc.recur(HG_row_m[i,], m[i], n[i], b_n[i], N, B)
    HG_row_n.either[i,] <- HG_row_n.calc(HG_row_m[i,], m[i], n[i], b_n[i], N, B)
  }
  
  expect_equal(HG_row_n.iter, HG_row_n.exp, info = "calculation using iteration", tolerance = 1e-4)
  expect_equal(HG_row_n.recur, HG_row_n.exp, info = "calculation using recursion", tolerance = 1e-4)
  expect_equal(HG_row_n.either, HG_row_n.exp, info = "calculation using either method", tolerance = 1e-4)  
})

mHG.statistic.calc.simple <- function(lambdas, n_max = length(lambdas)) {
  N <- length(lambdas)
  B <- sum(lambdas)
  b <- 0
  n <- 0
  
  mHG <- 1
  mHG.n <- n
  mHG.b <- b
  
  if (n_max > 0) {
    for (n in seq(n_max)) {
      b <- b + lambdas[n]
      HGT <- phyper((b-1), B, N - B, n, lower.tail = F)
      if (HGT < mHG) { # NOTE: Takes last, in order to be consistent
        mHG <- HGT
        mHG.n <- n
        mHG.b <- b
      }
    }
  }
  return(new("mHG.statistic.info", mHG = mHG, n = mHG.n, b = mHG.b))
}
generate.lambdas <- function(N, B) {
  lambdas <- numeric(N)
  lambdas[sample(N, B, replace = F)] <- 1
  return(lambdas)
}

test_that("mHG statistic calculation is consistent with R, n_max = N", {
  N.statistics <- 100
  N <- 100
  B <- 40
  
  for (i in seq(N.statistics)) {
    lambdas <- generate.lambdas(N, B)
    exp <- mHG.statistic.calc.simple(lambdas)
    calc <- mHG.statistic.calc(lambdas)
    expect_equal(exp, calc)
  }
})
test_that("mHG statistic calculation is consistent with R, B < n_max < N", {
  N.statistics <- 100
  N <- 100
  B <- 40
  n_max <- 50
  
  for (i in seq(N.statistics)) {
    lambdas <- generate.lambdas(N, B)
    exp <- mHG.statistic.calc.simple(lambdas, n_max = n_max)
    calc <- mHG.statistic.calc(lambdas, n_max = n_max)
    expect_equal(exp, calc)
  }
})
test_that("mHG statistic calculation is consistent with R, 0 < n_max < B", {
  N.statistics <- 100
  N <- 100
  B <- 40
  n_max <- 10
  
  for (i in seq(N.statistics)) {
    lambdas <- generate.lambdas(N, B)
    exp <- mHG.statistic.calc.simple(lambdas, n_max = n_max)
    calc <- mHG.statistic.calc(lambdas, n_max = n_max)
    expect_equal(exp, calc)
  }
})
test_that("mHG statistic calculation is consistent with R, n_max = 0", {
  N.statistics <- 100
  N <- 100
  B <- 40
  
  for (i in seq(N.statistics)) {
    lambdas <- generate.lambdas(N, B)
    exp <- mHG.statistic.calc.simple(lambdas, n_max = 0)
    expect_identical(exp, new("mHG.statistic.info", mHG = 1, n = 0, b = 0))
  }
})

test_that("mHG statistic calculation is consistent with R, n_max = 0", {
  N.statistics <- 100
  N <- 100
  B <- 40
  
  for (i in seq(N.statistics)) {
    lambdas <- generate.lambdas(N, B)
    exp <- mHG.statistic.calc.simple(lambdas, n_max = 0)
    expect_identical(exp, new("mHG.statistic.info", mHG = 1, n = 0, b = 0))
  }
})

test_that("mHG statistic calculation errors on invalid lambdas", {
  expect_error(mHG.statistic.calc(c()))
  expect_error(mHG.statistic.calc(c(1,1,0,2,2)))
})

test_that("mHG statistic calculation errors on invalid n_max", {
  N <- 100
  lambdas <- generate.lambdas(N = N, B = 40)
  expect_error(mHG.statistic.calc(lambdas, n_max = -1))
  expect_error(mHG.statistic.calc(lambdas, n_max = 0))
  expect_error(mHG.statistic.calc(lambdas, n_max = (N + 1)))
  expect_error(mHG.statistic.calc(lambdas, n_max = (N + 2)))
})


context("Pval Calculation")
no_R_separation_line <- function(B) {
  rep(0, B + 1)
}

expect_equal.pi_r <- function(pi_r, pi_r.exp, R_separation_line) {
  expect_equal(pi_r, pi_r.exp, info = paste("W (exp):", nrow(pi_r.exp) - 2,
                                            "B (exp):", ncol(pi_r.exp) - 2,
                                            "R separation line: ", paste(R_separation_line, collapse = ",")))
}

test_that("pi_r calculation, no R separation line, N = B", {
  for (i in 1:10) {
    N <- i
    B <- i
    R_separation_line <- no_R_separation_line(B)
    pi_r.exp <- matrix(0, nrow = N - B + 2,
                       ncol = B + 2, byrow = T)
    pi_r.exp[-1,-1] <- 1
    pi_r <- pi_r.calc(N, B, R_separation_line)
    expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
  }
})

test_that("pi_r calculation, no R separation line, predetermined scenarios", { 
  N <- 2
  B <- 1
  R_separation_line <- no_R_separation_line(B)
  pi_r.exp <- matrix(data = c(0,0,0,0,1,0.5,0,0.5,1), nrow = N - B + 2,
                     ncol = B + 2, byrow = T) # Follows from the definition of pi_r
  pi_r <- pi_r.calc(N, B, R_separation_line)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
  
  N <- 3
  B <- 2
  R_separation_line <- no_R_separation_line(B)
  pi_r.exp <- matrix(data = c(0,0,0,0,0,1,2/3,1/3,0,1/3,2/3,1), nrow = N - B + 2,
                     ncol = B + 2, byrow = T)
  pi_r <- pi_r.calc(N, B, R_separation_line)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
})

full_R_separation_line <- function(B, N) {
  rep(N - B, B + 1)
}
test_that("pi_r calculation, R separation line full", {
  N <- 10
  B <- 3
  R_separation_line <- full_R_separation_line(B, N)
  pi_r <- pi_r.calc(N, B, R_separation_line)
  pi_r.exp <- matrix(0, nrow = N - B + 2, ncol = B + 2)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
})

test_that("pi_r calculation, R separation line, predetermined scenarios", {
  N <- 3
  B <- 1
  R_separation_line <- c(0,2)
  pi_r.exp <- matrix(data = c(0,0,0,0,1,0,0,2/3,0,0,1/3,1/3), nrow = N - B + 2,
                     ncol = B + 2, byrow = T)
  pi_r <- pi_r.calc(N, B, R_separation_line)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
  
  N <- 3
  B <- 2
  R_separation_line <- c(0,1,1)
  pi_r.exp <- matrix(data = c(0,0,0,0,0,1,0,0,0,1/3,1/3,1/3), nrow = N - B + 2,
                     ncol = B + 2, byrow = T)
  pi_r <- pi_r.calc(N, B, R_separation_line)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
  
  N <- 4
  B <- 2
  R_separation_line <- c(0,1,1)
  pi_r.exp <- matrix(data = c(0,0,0,0,0,1,0,0,0,1/2,1/3,1/6,0,1/6,1/3,1/2), nrow = N - B + 2,
                     ncol = B + 2, byrow = T)
  pi_r <- pi_r.calc(N, B, R_separation_line)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
  
  R_separation_line <- c(0,1,2)
  pi_r.exp <- matrix(data = c(0,0,0,0,0,1,0,0,0,1/2,1/3,0,0,1/6,1/3,1/3), nrow = N - B + 2,
                     ncol = B + 2, byrow = T)
  pi_r <- pi_r.calc(N, B, R_separation_line)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
  
  R_separation_line <- c(0,2,2)
  pi_r.exp <- matrix(data = c(0,0,0,0,0,1,0,0,0,1/2,0,0,0,1/6,1/6,1/6), nrow = N - B + 2,
                     ncol = B + 2, byrow = T)
  pi_r <- pi_r.calc(N, B, R_separation_line)
  expect_equal.pi_r(pi_r, pi_r.exp, R_separation_line)
})

R_separation_line.test <- function(N, B, n_max = N, N.lines = 10) {
  W <- N - B
  for (i in seq(N.lines)) {
    HGT.mat <- matrix(nrow = W + 1, ncol = B + 1)
    for (w in 0:W) {
      for (b in 0:B) {
        HGT.mat[w + 1,b + 1] <- phyper((b-1), B, W, b + w, lower.tail = F)
      }
    }
    p <- HGT.mat[sample(length(HGT.mat), 1)]
    R_separation_line.exp <-  rep(W + 1, times = B + 1)
    
    for (b in 0:B) {
      p.obs <- -Inf
      w <- -1
      while ((p.obs <= p + EPSILON) && (w < W) && (w + b <= n_max)) {
        w <- w + 1
        p.obs <- HGT.mat[w + 1, b + 1]
      }
      if ((p.obs > p + EPSILON) || (w + b > n_max)) {
        R_separation_line.exp[b + 1] <- w
      } else {
        R_separation_line.exp[b + 1] <- W + 1        
      }
    }
    # Makes this a line
    highest_w_yet <- -1
    for (b in 0:B) {
      w <- R_separation_line.exp[b + 1]
      if (w > highest_w_yet) {
        highest_w_yet <- w
      } else {
        R_separation_line.exp[b + 1] <- highest_w_yet
      }
    }
    
    R_separation_line <- R_separation_line.calc(p, N, B, n_max = n_max)
    expect_equal(R_separation_line, R_separation_line.exp)
  }
}

test_that("R separation line consistent with R, n_max = N", {
  R_separation_line.test(20, 5)
})
test_that("R separation line, B < n_max < N", {
  R_separation_line.test(20, 5, n_max = 10)
})
test_that("R separation line, 0 < n_max < B", {
  R_separation_line.test(20, 5, n_max = 3)
})
test_that("R separation line, 0 < n_max < B", {
  R_separation_line.test(20, 5, n_max = 0)
})
test_that("R separation line, nmax = 0", {
  N <- 20
  B <- 5
  W <- N - B
  R_separation_line <- R_separation_line.calc(p = 0, 20, 5, n_max = 0)
  expect_equal(R_separation_line, numeric(length = B + 1))
  
  R_separation_line <- R_separation_line.calc(p = 1, 20, 5, n_max = 0)
  expect_equal(R_separation_line, numeric(length = B + 1) + 1)  
})

test_that("mHG.pval.calc, p < -EPSILON warning", {
  expect_warning(mHG.pval.calc(-(2 * EPSILON), 20, 5), 
                 "p-value calculation will be highly inaccurate due to an extremely small mHG statistic")  
})
test_that("mHG.pval.calc errors on invalid n_max", {
  N <- 10
  expect_error(mHG.pval.calc(p = 0.5, N = N, B = 5, n_max = -1))
  expect_error(mHG.pval.calc(p = 0.5, N = N, B = 5, n_max = 0))
  expect_error(mHG.pval.calc(p = 0.5, N = N, B = 5, n_max = N + 1)) 
  expect_error(mHG.pval.calc(p = 0.5, N = N, B = 5, n_max = N + 2))
})
test_that("mHG.pval.calc errors when B > N", {
  N <- 10
  expect_error(mHG.pval.calc(p = 0.5, N = N, B = N + 1))
  expect_error(mHG.pval.calc(p = 0.5, N = N, B = N + 2))
})
test_that("mHG.pval.calc errors when p > 1", {
  N <- 10
  expect_error(mHG.pval.calc(p = 1.1, N = N, B = N + 1))
})

mHG.pval.calc.test <- function(N, B,n_max = N, mHG = c(), p.exp = c()) {
  mHG <- c(-EPSILON, mHG, 1)
  p.exp <- c(0, p.exp, 1)
  p <- sapply(mHG, function(x) mHG.pval.calc(x, N, B, n_max = n_max))
  expect_equal(p, p.exp)
}

# Code used to generate such scenarios:
#   N <- ...
#   B <- ...
#   combn.res <- combn(N, B)
#   mHG <- numeric(ncol(combn.res))
#   for (i in 1:ncol(combn.res)) {
#     lambdas <- numeric(N)
#     lambdas[combn.res[,i]] <- 1
#     mHG[i] <- mHG.statistic.calc(lambdas)@mHG
#   }
#   table(mHG)
test_that("mHG.pval.calc, n_max = N, predetermined scenarios", {
  mHG.pval.calc.test(N = 3, B = 3)
  mHG.pval.calc.test(N = 4, B = 2, mHG = c(1/6,1/2,5/6), p.exp = c(1/6,2/3,5/6))
  mHG.pval.calc.test(N = 4, B = 3, mHG = c(1/4,1/2,3/4), p.exp = c(1/4,1/2,3/4))
  mHG.pval.calc.test(N = 5, B = 3, mHG = c(1/10,3/10,2/5,3/5,7/10,9/10), p.exp = c(1/10,3/10,1/2,7/10,8/10,9/10))
  mHG.pval.calc.test(N = 6, B = 2, mHG = c(1/15,1/5,1/3,2/5,3/5,2/3,4/5,14/15), p.exp = c(1/15,1/5,2/5,8/15,2/3,4/5,13/15,14/15))
  mHG.pval.calc.test(N = 6, B = 3, mHG = c(1/20,1/5,1/2,4/5,19/20), p.exp = c(1/20,3/10,3/4,9/10,19/20))
})
test_that("mHG.pval.calc, n_max < N, predetermined scenarios", {
  mHG.pval.calc.test(N = 6, B = 2, n_max = 1, mHG = c(1/3), p.exp = c(1/3))
  mHG.pval.calc.test(N = 6, B = 2, n_max = 2, mHG = c(1/15,1/3,3/5), p.exp = c(1/15,1/3,3/5))
  mHG.pval.calc.test(N = 6, B = 2, n_max = 3, mHG = c(1/15,1/5,1/3,3/5,4/5), p.exp = c(1/15,1/5,2/5,3/5,4/5))
  mHG.pval.calc.test(N = 6, B = 2, n_max = 4, mHG = c(1/15,1/5,1/3,2/5,3/5,4/5,14/15), p.exp = c(1/15,1/5,2/5,8/15,2/3,4/5,14/15))
  mHG.pval.calc.test(N = 6, B = 2, n_max = 5, mHG = c(1/15,1/5,1/3,2/5,3/5,2/3,4/5,14/15), p.exp = c(1/15,1/5,2/5,8/15,2/3,4/5,13/15,14/15))
  mHG.pval.calc.test(N = 6, B = 2, n_max = 6, mHG = c(1/15,1/5,1/3,2/5,3/5,2/3,4/5,14/15), p.exp = c(1/15,1/5,2/5,8/15,2/3,4/5,13/15,14/15))
})

context("mHG Test")

mHG.test.test <- function(lambdas, n_max = length(lambdas), mHG.exp, n.exp, p.exp) {
  htest <- mHG.test(lambdas = lambdas, n_max = n_max)
  expect_is(htest, "htest")
  expect_equal(unname(htest$statistic), mHG.exp)
  expect_equal(unname(htest$parameters["N"]), length(lambdas))
  expect_equal(unname(htest$parameters["B"]), sum(lambdas))
  expect_equal(unname(htest$parameters["n_max"]), n_max)
  expect_equal(htest$n, n.exp)
  if (n.exp > 0) {
    expect_equal(htest$b, sum(lambdas[1:n.exp]))
  } else {
    expect_equal(htest$b, 0)
  }
  expect_equal(htest$p.value, p.exp)
}

test_that("mHG.test work, n_max = N, predetermined scenarios", {
  mHG.test.test(lambdas = c(0,0), mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,1), mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(1,0), mHG.exp = 1/2, n.exp = 1, p.exp = 1/2)
  mHG.test.test(lambdas = c(1,1), mHG.exp = 1, n.exp = 0, p.exp = 1)
  
  mHG.test.test(lambdas = c(0,0,0), mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,0,1), mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,1,0), mHG.exp = 2/3, n.exp = 2, p.exp = 2/3)
  mHG.test.test(lambdas = c(0,1,1), mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(1,0,0), mHG.exp = 1/3, n.exp = 1, p.exp = 1/3)
  mHG.test.test(lambdas = c(1,0,1), mHG.exp = 2/3, n.exp = 1, p.exp = 2/3)
  mHG.test.test(lambdas = c(1,1,0), mHG.exp = 1/3, n.exp = 2, p.exp = 1/3)
  mHG.test.test(lambdas = c(1,1,1), mHG.exp = 1, n.exp = 0, p.exp = 1)  

})

test_that("mHG.test work, n_max < N, predetermined scenarios", {
  mHG.test.test(lambdas = c(0,0,0), n_max = 1, mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,0,1), n_max = 1, mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,1,0), n_max = 1, mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,1,1), n_max = 1, mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(1,0,0), n_max = 1, mHG.exp = 1/3, n.exp = 1, p.exp = 1/3)
  mHG.test.test(lambdas = c(1,0,1), n_max = 1, mHG.exp = 2/3, n.exp = 1, p.exp = 2/3)
  mHG.test.test(lambdas = c(1,1,0), n_max = 1, mHG.exp = 2/3, n.exp = 1, p.exp = 2/3)
  mHG.test.test(lambdas = c(1,1,1), n_max = 1, mHG.exp = 1, n.exp = 0, p.exp = 1)  
  
  mHG.test.test(lambdas = c(0,0,0), n_max = 2, mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,0,1), n_max = 2, mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(0,1,0), n_max = 2, mHG.exp = 2/3, n.exp = 2, p.exp = 2/3)
  mHG.test.test(lambdas = c(0,1,1), n_max = 2, mHG.exp = 1, n.exp = 0, p.exp = 1)
  mHG.test.test(lambdas = c(1,0,0), n_max = 2, mHG.exp = 1/3, n.exp = 1, p.exp = 1/3)
  mHG.test.test(lambdas = c(1,0,1), n_max = 2, mHG.exp = 2/3, n.exp = 1, p.exp = 2/3)
  mHG.test.test(lambdas = c(1,1,0), n_max = 2, mHG.exp = 1/3, n.exp = 2, p.exp = 1/3)
  mHG.test.test(lambdas = c(1,1,1), n_max = 2, mHG.exp = 1, n.exp = 0, p.exp = 1)  
})