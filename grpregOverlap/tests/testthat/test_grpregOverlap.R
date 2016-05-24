library(testthat)
library(grpregOverlap)

context("Testing grpregOverlap()")

test_that("Non-overlapping fit againt grpreg:", {
  data(birthwt.grpreg)
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  ## no overlap, should be the same as from 'grpreg'
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(9), c(10, 11), 
                c(12), c(13), c(14, 15, 16))
  group2 <- rep(1:length(group), sapply(group, length))
  
  ## linear regression
  y <- birthwt.grpreg$bwt
  invisible(capture.output({
    fit <- grpregOverlap(X, y, group, penalty = 'grLasso')
  }))
  # x == x.latent, test expandX()
  expect_equal(all(X == fit$X.latent), TRUE) 
  # beta = beta.latent, test gamma2beta()
  # Pass test ONLY IF variable indices in ascending order, otherwise, 
  # variables are reordered in beta.latent according to group order.
  expect_equal(all(fit$beta == fit$beta.latent), TRUE)
  
  # equivalent to fit 'grpreg'
  fit2 <- grpreg(X, y, group2, penalty = 'grLasso')
  expect_identical(fit$beta, fit2$beta)
  expect_equal(all(fit2$beta == fit$beta.latent), TRUE)
  
  ## logistic regression
  y <- birthwt.grpreg$low
  invisible(capture.output({
    fit <- grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')
  }))
  # x == x.latent, test expandX()
  expect_equal(all(X == fit$X.latent), TRUE) 
  # beta = beta.latent, test gamma2beta()
  # Pass test ONLY IF variable indices in ascending order, otherwise, 
  # variables are reordered in beta.latent according to group order.
  expect_equal(all(fit$beta == fit$beta.latent), TRUE)
  
  # equivalent to fit 'grpreg'
  fit2 <- grpreg(X, y, group2, penalty = 'grLasso', family = 'binomial')
  expect_identical(fit$beta, fit2$beta)
  expect_equal(all(fit2$beta == fit$beta.latent), TRUE)
  
  ## logistic regression
  y <- birthwt.grpreg$low
  invisible(capture.output({
    fit <- grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')
  }))
  # x == x.latent, test expandX()
  expect_equal(all(X == fit$X.latent), TRUE) 
  # beta = beta.latent, test gamma2beta()
  # Pass test ONLY IF variable indices in ascending order, otherwise, 
  # variables are reordered in beta.latent according to group order.
  expect_equal(all(fit$beta == fit$beta.latent), TRUE)
  
  # equivalent to fit 'grpreg'
  fit2 <- grpreg(X, y, group2, penalty = 'grLasso', family = 'binomial')
  expect_identical(fit$beta, fit2$beta)
  expect_equal(all(fit2$beta == fit$beta.latent), TRUE)
})


test_that("predict, coef, select, cv, against grpreg: ", {
  data(birthwt.grpreg)
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  ## no overlap, should be the same as from 'grpreg'
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(9), c(10, 11), 
                c(12), c(13), c(14, 15, 16))
  group2 <- rep(1:length(group), sapply(group, length))
  ## logistic regression
  y <- birthwt.grpreg$low
  invisible(capture.output({
    fit <- grpregOverlap(X, y, group, penalty = 'grMCP', family = 'binomial')
  }))
  
  ## test predict, coef againt grpreg
  fit2 <- grpreg(X, y, group2, penalty = 'grMCP', family = 'binomial')
  
  expect_equal(all(coef(fit, lambda=.001) == coef(fit2, lambda=.001)), TRUE)
  expect_equal(all(predict(fit, X, type="link", lambda=.001) == 
                     predict(fit2, X, type="link", lambda=.001)
    ), TRUE)
  expect_equal(all(predict(fit, X, type="response", lambda=.001) ==
                     predict(fit2, X, type="response", lambda=.001)
    ), TRUE)
  expect_equal(all(
    ), TRUE)
  
  expect_equal(all(predict(fit, X, type="class", lambda=.001) ==
                     predict(fit2, X, type="class", lambda=.001)
    ), TRUE)
  expect_equal(all(predict(fit, type="vars", lambda=.07) ==
                     predict(fit, type='vars', latent=T, lambda=.07)
    ), TRUE)
  expect_equal(all(predict(fit, type="nvars", lambda=.07) ==
                     predict(fit, type='nvars', latent=T, lambda=.07)
  ), TRUE)
  expect_equal(all(predict(fit, type="coefficients", lambda=.07) ==
                     predict(fit, type="coefficients", latent = T, lambda=.07)
    ), TRUE)

  expect_equal(all(predict(fit, type="vars", lambda=.07) ==
                     predict(fit2, type="vars", lambda=.07)
    ), TRUE)
  
  expect_equal(all(predict(fit, type="groups", lambda=.07) ==
                     predict(fit2, type="groups", lambda=.07)
    ), TRUE)
  
  invisible(capture.output({
    expect_equal(all(predict(fit, type="norm", lambda=.07) ==
                       predict(fit2, type="norm", lambda=.07)
    ), TRUE)
    }))
  
  invisible(capture.output({
    cvfit <- cv.grpregOverlap(X, y, group, family="binomial", penalty="grMCP", 
                               seed = 1234)
  }))
  
  cvfit2 <- cv.grpreg(X, y, group2, family="binomial", penalty="grMCP", 
                      seed = 1234)
  expect_equal(all(coef(cvfit) == coef(cvfit2)), TRUE)
  expect_equal(all(predict(cvfit, X) == predict(cvfit2, X)), TRUE)
  expect_equal(all(predict(cvfit, X, type="response") ==
                     predict(cvfit2, X, type="response")
    ), TRUE)
  expect_equal(all(predict(cvfit, type="groups") == 
                     predict(cvfit2, type="groups")
    ), TRUE)
  
  ## test select
  y <- birthwt.grpreg$bwt
  invisible(capture.output({
    fit <- grpregOverlap(X, y, group, penalty="grLasso")
  }))
  fit2 <- grpreg(X, y, group2, penalty="grLasso")
  expect_equal(all(select(fit)$lambda == select(fit2)$lambda), TRUE)
  expect_equal(all(select(fit)$beta == select(fit2)$beta), TRUE)
  suppressWarnings(
    expect_equal(all(select(fit,crit="AIC",df="active")$beta.latent ==
                     select(fit2,crit="AIC",df="active")$beta
     ),TRUE)  
  )  
})
# 
# test_that("Overlapping fit: ", {
#   
#   ## linear regression, a simulation demo.
#   set.seed(123)
#   group <- list(gr1 = c(1, 2, 3),
#                 gr2 = c(1, 4),
#                 gr3 = c(2, 4, 5),
#                 gr4 = c(3, 5),
#                 gr5 = c(6))
#   
#   beta.latent.T <- c(5, 5, 5, 0, 0, 0, 0, 0, 5, 5, 0) # true latent coefficients.
#   # beta.T <- c(2, 3, 7, 0, 5, 0), true variables: 1, 2, 3, 5; true groups: 1, 4.
#   X <- matrix(rnorm(n = 6*100), ncol = 6)
#   incid.mat <- incidenceMatrix(X, group) # group membership incidence matrix
#   over.mat <- Matrix(incid.mat %*% t(incid.mat)) # overlap matrix
# #   grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector  
#   X.latent <- expandX(X, group)
#   y <- X.latent %*% beta.latent.T + rnorm(100)
#   
#   fit <- grpregOverlap(X, y, group, penalty = 'grLasso')
# #   fit <- grpregOverlap(X, y, group, penalty = 'grMCP')
# #   fit <- grpregOverlap(X, y, group, penalty = 'grSCAD')
#   coef(fit, latent = TRUE) # compare to beta.latent.T
#   plot(fit, latent = TRUE) 
#   coef(fit, latent = FALSE) # compare to beta.T
#   plot(fit, latent = FALSE)
#   predict(fit, type = 'vars', latent = T, lambda = 4)
#   predict(fit, type = 'vars', latent = F, lambda = 4)
#   select(fit, "BIC")
#   select(fit, "AIC")
# 
#   cvfit <- cv.grpregOverlap(X, y, group, penalty = 'grMCP')
#   plot(cvfit)
#   par(mfrow=c(2,2))
#   plot(cvfit, type="all")
#   coef(cvfit)
#   summary(cvfit)
# 
#   
#   ## logistic regression, real data, pathway selection
#   data(pathway.dat)
#   X <- pathway.dat$expression
#   group <- pathway.dat$pathways
#   y <- pathway.dat$mutation
#   fit <- grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')
# #   fit <- grpregOverlap(X, y, group, penalty = 'grMCP', family = 'binomial')
# #   fit <- grpregOverlap(X, y, group, penalty = 'grSCAD', family = 'binomial')
# #   fit <- grpregOverlap(X, y, group, penalty = 'gel', family = 'binomial')
# #   fit <- grpregOverlap(X, y, group, penalty = 'cMCP', family = 'binomial')
# #   print(object.size(cvfit), units = 'Mb')
#   plot(fit)
#   plot(fit, latent = FALSE)
#   plot(fit, norm = TRUE)
#   plot(fit, norm = T)
# 
#   predict(fit, type = 'ngroups', lambda = 0.01)
#   predict(fit, type = 'nvars', lambda = 0.01)
#   predict(fit, type = 'vars', latent = TRUE, lambda = 0.01)
#   predict(fit, type = 'groups', latent = TRUE, lambda = 0.01) # A note printed.
#   predict(fit, X, type="class", lambda=.01)
#   predict(fit, X, type = "coefficients", lambda = 0.01)
#   predict(fit, type="norm", lambda=.01)
#   
#   select(fit)
#   select(fit,crit="AIC",df="active")
# 
#   cvfit <- cv.grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')
#   coef(cvfit)
#   predict(cvfit, X, type='response')
#   predict(cvfit, X, type = 'class')
#   plot(cvfit)
#   plot(cvfit, type = 'all')
#   summary(cvfit)
# 
# })
