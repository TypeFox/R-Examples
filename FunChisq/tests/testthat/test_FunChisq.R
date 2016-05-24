library(testthat)
library(FunChisq)
context("Testing the FunChisq package")


test_that("Testing the exact functional test", {

  exact.functional.test <- ExactFunctionalTest

  x1 <- matrix(c(12, 26, 18, 0, 8, 12), nrow=2, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x1), 8), 0.042556227)
  expect_equal(signif(exact.functional.test(t(x1)), 8), 0.027271581)

  x2 <- matrix(c(0,0,0,0,0,0,0,0,0), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x2), 1)

  x3 <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  expect_equal(signif(exact.functional.test(x3), 8), 0.002997003)
  expect_equal(signif(exact.functional.test(t(x3)), 8), 0.0065490065)

  if(0) { # This test case causes hang on windows. To be fixed.
    x4 <- matrix(rep(10,25), nrow=5)
    expect_equal(exact.functional.test(x4), 1)
  }

  x5 <- matrix(c(4,0,0,0,4,0,0,0,4), nrow=3, byrow = TRUE)
  expect_equal(signif(exact.functional.test(x5), 8), 0.00017316017)

  x6 <- matrix(c(2,0,0,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x6), 8),
               signif(stats::fisher.test(x6)$p.value, 8))

  x7 <- matrix(c(2,2,2,2), nrow=2, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x7), 8),
               signif(stats::fisher.test(x7)$p.value, 8))

  x8 <- matrix(c(0,10,15,20,5,0,25,0,0), nrow=3, byrow = TRUE)
  expect_equivalent(signif(exact.functional.test(x8), 8),
               signif(fun.chisq.test(x8)$p.value, 8))

  x9 <- matrix(c(1,1,1,1,1,1,1,1,1), nrow=3, byrow = TRUE)
  expect_equal(exact.functional.test(x9), 1)
})


test_that("Testing the functional chi-square test", {

  ex <- list()

  ex[[1]] <- list(
    x = matrix(c(20,0,20,0,20,0,5,0,5), 3),
    method = "default",
    stat.truth = 72,
    par.truth = 4,
    estimate.cond.fun.index = 0.77459667,
    estimate.fun.index = 0.71713717,
    pval.truth = 8.5822345e-15,
    digits = 8
  )

  ex[[2]] <- list(
    x = ex[[1]]$x,
    method = "nfchisq",
    stat.truth = 24.0416306,
    par.truth = 4,
    estimate.cond.fun.index = 0.77459667,
    estimate.fun.index = 0.71713717,
    pval.truth = 5.1061401e-128,
    digits = 8
  )

  ex[[3]] <- list(
    x = t(ex[[1]]$x),
    method = "fchisq",
    stat.truth = 64.285714,
    par.truth = 4,
    estimate.cond.fun.index = 0.67936622,
    estimate.fun.index = 0.67763093,
    pval.truth = 3.6385174e-13,
    digits = 8
  )

  ex[[4]] <- list(
    x = ex[[3]]$x,
    method = "nfchisq",
    stat.truth = 21.314219,
    par.truth = 4,
    estimate.cond.fun.index = 0.67936622,
    estimate.fun.index = 0.67763093,
    pval.truth = 4.1897164e-101,
    digits = 8
  )

  ex[[5]] <- list(
    x = matrix(c(5, 0, 0, 0, 5, 0, 0, 0, 5), nrow=3),
    method = "fchisq",
    stat.truth = 30,
    par.truth = 4,
    estimate.cond.fun.index = 1,
    estimate.fun.index = 1,
    pval.truth = 4.8944371e-06,
    digits = 8
  )

  ex[[6]] <- list(
    x = matrix(c(5, 0, 5, 0, 0, 5, 0, 5), nrow=4),
    method = "fchisq",
    stat.truth = 20,
    par.truth = 3,
    estimate.cond.fun.index = 1,
    estimate.fun.index = 1,
    pval.truth = 0.00016974244,
    digits = 8
  )

  for(i in seq_along(ex)) {
    within(ex[[i]],
           {
             h <- fun.chisq.test(x, method=method)
             expect_equivalent(h$statistic, stat.truth)
             expect_equivalent(h$parameter, par.truth)
             expect_equivalent(signif(h$estimate, digits=digits), estimate.fun.index)
             expect_equivalent(signif(h$p.value, digits=digits), pval.truth)

             h <- fun.chisq.test(x, method=method, index.kind="conditional")
             expect_equivalent(signif(h$estimate, digits=digits), estimate.cond.fun.index)
           }
    )
  }
})


test_that("Testing the comparative functional chi-square test", {
  x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  y <- t(x)
  z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
  data <- list(x,y,z)
  expect_equivalent(signif(cp.fun.chisq.test(data)$p.value, 8), 0.00018762119)
  expect_equivalent(signif(cp.fun.chisq.test(data, method="nfchisq")$p.value, 8),
               1.0052639e-07)
})


test_that("Testing the comparative chi-square test", {

  x <- list()

  x[[1]] <- matrix(c(0,0,0,
                     0,0,0,
                     0,0,0), nrow=3)
  x[[2]] <- x[[1]]
  x[[3]] <- x[[1]]

  h <- cp.chisq.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 0)

  h <- cp.chisq.test(x, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 1)

  x <- list()

  x[[1]] <- matrix(c(4,0,0,
                     0,4,0,
                     0,0,4), nrow=3)
  x[[2]] <- x[[1]]
  x[[3]] <- x[[1]]
  h <- cp.chisq.test(x)
  expect_equivalent(signif(h$p.value, 8), 1)
  expect_equivalent(signif(h$statistic, 8), 0)
  expect_equivalent(h$parameter, 8)

  h <- cp.chisq.test(x, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 0.97724987)

  x <- list()

  x[[1]] <- matrix(c(4,0,0,
                     0,4,0,
                     0,0,4), nrow=3)

  x[[2]] <- matrix(c(0,4,4,
                     4,0,4,
                     4,4,0), nrow=3)

  h <- cp.chisq.test(x)
  expect_equivalent(signif(h$p.value, 8), 2.8936962e-07)
  expect_equivalent(signif(h$statistic, 8), 36)
  expect_equivalent(h$parameter, 4)

  h <- cp.chisq.test(x, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 0)

  x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
  y <- t(x)
  z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
  data <- list(x,y,z)
  h <- cp.chisq.test(data)
  expect_equivalent(signif(h$p.value, 8), 1.3542453e-06)
  expect_equivalent(signif(h$statistic, 8), 42)
  expect_equivalent(h$parameter, 8)

  h <- cp.chisq.test(data, method="nchisq")
  expect_equivalent(signif(h$p.value, 8), 9.4795348e-18)
})


test_that("Testing the interaction analysis", {
  x <- matrix(c(0,0,1,0,1,
                1,0,2,1,0,
                2,2,0,0,0,
                1,2,1,1,2,
                1,0,2,1,2),
              nrow = 5, ncol = 5, byrow = TRUE)
  #x <- data.frame(x)
  independent.variables <- list(c(1),c(1),c(1),
                              c(2),c(2),c(2),
                              c(1,2), c(2,3),
                              c(3,4), c(4,5),
                              c(1,2,3), c(1,2,3,4))
  dependent.variable <- c(3,4,5,
                            3,4,5,
                            3,4,
                            5,1,
                            4,5)
  #variable.names <- c("variable1", "variable2", "variable3", "variable4", "variable5")

  #randon variables.
  #But quantized expressions should be not all the same.
  #
  # x <- t(sapply(c(1:10), function(y){
  #   sample(x = c(1:3), size = 6, replace = TRUE)
  #   return()
  # }))
  # #x <- data.frame(x)
  # dependent.variables <- lapply(c(1,1,2,2,3,3), function(y){
  #   sample(c(1:nrow(x)), size = y)
  # })
  # independent.variable <- sample(c(1:nrow(x)), size = 6, replace = TRUE)

  index.kind <- "unconditional"

  output.C <- test.interactions(x = x,
                                list.ind.vars = independent.variables,
                                dep.var = dependent.variable,
                                index.kind = index.kind)

  output.C <- output.C[, c(3:ncol(output.C))]

  output.R <- NULL
  for(i in c(1:length(dependent.variable))){
    P <- x[independent.variables[[i]],]
    C <- x[dependent.variable[[i]],]

    output.tmp <- NULL
    if(is.vector(P)){
      t <- table(P, C)
      fc <- fun.chisq.test(t, index.kind = "unconditional")
      output.tmp <- c(fc$p.value, fc$statistic, fc$estimate)
    }else{
      A <- c(lapply(c(1:nrow(P)), function(y){return(unlist(P[y,]))}), list(unlist(C)))
      A <- ftable(A, row.vars = 1:(length(A)-1), col.vars = length(A))

      A <- A[!sapply(c(1:nrow(A)), function(y){all(A[y,]==0)}),] # remove all 0 rows
      A <- A[,!sapply(c(1:ncol(A)), function(y){all(A[,y]==0)})] # remove all 0 columns

      fc <- fun.chisq.test(A, index.kind = "unconditional")
      output.tmp <- c(fc$p.value, fc$statistic, fc$estimate)
    }

    if(is.null(output.R)){
      output.R <- output.tmp
    }else{
      output.R <- rbind(output.R, output.tmp)
    }
  }
  colnames(output.R) <- c("p.value", "statistic", "estimate")
  rownames(output.R) <- c(1:nrow(output.R))
  output.R <- data.frame(output.R)
  #output.R <- data.frame(cbind(Parent=N[I[,1]], Child=N[I[,2]], output.R))

  output.diff <- output.C - output.R
  for(i in c(1:nrow(output.diff))){
    for(j in c(1:ncol(output.diff))){
      expect_equivalent(as.numeric(output.diff[i,j]), 0)
    }
  }
  # apply(output.diff, c(1,2), function(x){
  #   expect_equivalent(as.numeric(x), 0)
  #   return(NULL)
  # })
})
