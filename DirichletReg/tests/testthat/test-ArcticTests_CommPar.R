#TO DO:
#anova(resC1)
#residuals(resC1)
#update(resC1)



cat("================================================================================\n")
cat("=== Arctic Lake Data Tests =====================================================\n")
cat("================================================================================\n")

context("  AL: Original Data\n   ")

test_that("Arctic Lake - Original Data Structure", {
  expect_true(exists("ArcticLake"))
  expect_identical(dim(ArcticLake), c(39L, 4L))
  expect_identical(names(ArcticLake), c("sand", "silt", "clay", "depth"))
  expect_true(all(unlist(lapply(ArcticLake, function(colElement){ class(colElement) == "numeric" }))))
})

context("  AL: Transformation\n   ")

AL <- ArcticLake[, 4, drop=FALSE]

test_that("Arctic Lake - Data Transformation", {
  expect_identical(dim(AL), c(39L, 1L))
  expect_identical(class(AL), "data.frame")
  expect_warning(AL$Y <<- DR_data(ArcticLake[, 1:3]), ".*normalization forced.*")
  expect_equal(unname(rowSums(AL$Y)), rep(1.0, 39L))
})

cat("\n  --- Checks: Common Model -----------------------------------------------------")

#
#
#
#
#

context("  AL: Common - Null Model ( Y ~ 1 )\n   ")

resC1 <- DirichReg(Y ~ 1, AL)

load("testdata/resC1.RData")

test_that("Model Estimation", {
  expect_equal(resC1_mathematica$MLE, resC1$logLik)
  expect_equal(resC1_mathematica$DEV, -2.0*resC1$logLik)
  expect_equal(resC1_mathematica$COEFS, unname(resC1$coefficients), check.attributes = FALSE)
  expect_equal(resC1_mathematica$SE, unname(resC1$se), check.attributes = FALSE)
  expect_equal(resC1_mathematica$Z, unname(resC1$coefficients / resC1$se), check.attributes = FALSE)
  expect_equal(resC1_mathematica$P, 2*pnorm(-abs(unname(resC1$coefficients / resC1$se))), check.attributes = FALSE)
  expect_equal(resC1_mathematica$HESSIAN, unname(resC1$hessian), check.attributes = FALSE)
  expect_equal(resC1_mathematica$VCOV, unname(resC1$vcov), check.attributes = FALSE)
})

test_that("Methods", {
  expect_equal(resC1_mathematica$NOBS, nobs(resC1), check.attributes = FALSE)
  expect_equal(resC1_mathematica$MLE, unclass(logLik(resC1)), check.attributes = FALSE)
  expect_equal(resC1_mathematica$NPAR, attributes(logLik(resC1))$df, check.attributes = FALSE)
  expect_equal(resC1_mathematica$AIC, AIC(resC1), check.attributes = FALSE)
  expect_equal(resC1_mathematica$BIC, BIC(resC1), check.attributes = FALSE)
  expect_equal(resC1_mathematica$COEFS, unlist(coef(resC1)), check.attributes = FALSE)
  expect_equal(resC1_mathematica$VCOV, vcov(resC1), check.attributes = FALSE)

  expect_equal(resC1_mathematica$PREDICT$ALPHA, unname(fitted(resC1, alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
  expect_equal(resC1_mathematica$PREDICT$PHI,   unname(fitted(resC1, alpha=F, phi=T, mu=F))[1], check.attributes = FALSE)
  expect_equal(resC1_mathematica$PREDICT$MU,    unname(fitted(resC1, alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  expect_equal(resC1_mathematica$PREDICT$ALPHA, unname(predict(resC1, data.frame("depth"=0), alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
  expect_equal(resC1_mathematica$PREDICT$PHI,   unname(predict(resC1, data.frame("depth"=0), alpha=F, phi=T, mu=F)), check.attributes = FALSE)
  expect_equal(resC1_mathematica$PREDICT$MU,    unname(predict(resC1, data.frame("depth"=0), alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  conf_ints <- confint(resC1, level = c(.99, .95))
  conf_ints <- lapply(1:3, function(listelement){ sort(unlist(lapply(c(list(conf_ints$coefficients), conf_ints$ci), `[[`, listelement))) })
  conf_ints <- unname(t(matrix(unlist(conf_ints), 5)))
  expect_equal(resC1_mathematica$CONFINT, conf_ints)
})

#
#
#
#
#

context("  AL: Common - Linear ( Y ~ depth )\n   ")

resC2 <- DirichReg(Y ~ depth, AL)

load("testdata/resC2.RData")

test_that("Model Estimation", {
  expect_equal(resC2_mathematica$MLE, resC2$logLik)
  expect_equal(resC2_mathematica$DEV, -2.0*resC2$logLik)
  expect_equal(resC2_mathematica$COEFS, unname(resC2$coefficients), check.attributes = FALSE)
  expect_equal(resC2_mathematica$SE, unname(resC2$se), check.attributes = FALSE)
  expect_equal(resC2_mathematica$Z, unname(resC2$coefficients / resC2$se), check.attributes = FALSE)
  expect_equal(resC2_mathematica$P, 2*pnorm(-abs(unname(resC2$coefficients / resC2$se))), check.attributes = FALSE)
  expect_equal(resC2_mathematica$HESSIAN, unname(resC2$hessian), check.attributes = FALSE)
  expect_equal(resC2_mathematica$VCOV, unname(resC2$vcov), check.attributes = FALSE)
})

test_that("Methods", {
  expect_equal(resC2_mathematica$NOBS, nobs(resC2), check.attributes = FALSE)
  expect_equal(resC2_mathematica$MLE, unclass(logLik(resC2)), check.attributes = FALSE)
  expect_equal(resC2_mathematica$NPAR, attributes(logLik(resC2))$df, check.attributes = FALSE)
  expect_equal(resC2_mathematica$AIC, AIC(resC2), check.attributes = FALSE)
  expect_equal(resC2_mathematica$BIC, BIC(resC2), check.attributes = FALSE)
  expect_equal(resC2_mathematica$COEFS, unlist(coef(resC2)), check.attributes = FALSE)
  expect_equal(resC2_mathematica$VCOV, vcov(resC2), check.attributes = FALSE)

#  expect_equal(resC2_mathematica$PREDICT$ALPHA, unname(fitted(resC2, alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
#  expect_equal(resC2_mathematica$PREDICT$PHI,   unname(fitted(resC2, alpha=F, phi=T, mu=F))[1], check.attributes = FALSE)
#  expect_equal(resC2_mathematica$PREDICT$MU,    unname(fitted(resC2, alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  expect_equal(resC2_mathematica$PREDICT$ALPHA, unname(predict(resC2, data.frame("depth"=0:150), alpha=T, phi=F, mu=F)), check.attributes = FALSE)
  expect_equal(resC2_mathematica$PREDICT$PHI,   unname(predict(resC2, data.frame("depth"=0:150), alpha=F, phi=T, mu=F)), check.attributes = FALSE)
  expect_equal(resC2_mathematica$PREDICT$MU,    unname(predict(resC2, data.frame("depth"=0:150), alpha=F, phi=F, mu=T)), check.attributes = FALSE)

  # BEWARE! UGLY CODE AHEAD!
  conf_ints <- confint(resC2, level = c(.99, .95))
  conf_mat <- matrix(NA, 6, 4); rowz <- rep(1:2, 3); listz <- rep(1:3, each=2)
  for(i in 1:6){
    conf_mat[i,] <- c(conf_ints$ci[[2L]][[listz[i]]][rowz[i], 1L], conf_ints$ci[[1L]][[listz[i]]][rowz[i], ], conf_ints$ci[[2L]][[listz[i]]][rowz[i], 2L])
  }
  conf_ints <- unname(cbind(conf_mat[,1:2], unlist(conf_ints$coefficients), conf_mat[,3:4]))
  expect_equal(resC2_mathematica$CONFINT, conf_ints)
})

#
#
#
#
#

context("  AL: Common - Linear + Quadratic ( Y ~ depth + I(depth^2) )\n   ")

resC3 <- DirichReg(Y ~ depth + I(depth^2), AL)

load("testdata/resC3.RData")

test_that("Model Estimation", {
  expect_equal(resC3_mathematica$MLE, resC3$logLik)
  expect_equal(resC3_mathematica$DEV, -2.0*resC3$logLik)
  expect_equal(resC3_mathematica$COEFS, unname(resC3$coefficients), check.attributes = FALSE)
  expect_equal(resC3_mathematica$SE, unname(resC3$se), check.attributes = FALSE)
  expect_equal(resC3_mathematica$Z, unname(resC3$coefficients / resC3$se), check.attributes = FALSE)
  expect_equal(resC3_mathematica$P, 2*pnorm(-abs(unname(resC3$coefficients / resC3$se))), check.attributes = FALSE)
  expect_equal(resC3_mathematica$HESSIAN, unname(resC3$hessian), check.attributes = FALSE)
  expect_equal(resC3_mathematica$VCOV, unname(resC3$vcov), check.attributes = FALSE)
})

test_that("Methods", {
  expect_equal(resC3_mathematica$NOBS, nobs(resC3), check.attributes = FALSE)
  expect_equal(resC3_mathematica$MLE, unclass(logLik(resC3)), check.attributes = FALSE)
  expect_equal(resC3_mathematica$NPAR, attributes(logLik(resC3))$df, check.attributes = FALSE)
  expect_equal(resC3_mathematica$AIC, AIC(resC3), check.attributes = FALSE)
  expect_equal(resC3_mathematica$BIC, BIC(resC3), check.attributes = FALSE)
  expect_equal(resC3_mathematica$COEFS, unlist(coef(resC3)), check.attributes = FALSE)
  expect_equal(resC3_mathematica$VCOV, vcov(resC3), check.attributes = FALSE)

#  expect_equal(resC3_mathematica$PREDICT$ALPHA, unname(fitted(resC3, alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
#  expect_equal(resC3_mathematica$PREDICT$PHI,   unname(fitted(resC3, alpha=F, phi=T, mu=F))[1], check.attributes = FALSE)
#  expect_equal(resC3_mathematica$PREDICT$MU,    unname(fitted(resC3, alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  expect_equal(resC3_mathematica$PREDICT$ALPHA, unname(predict(resC3, data.frame("depth"=0:150), alpha=T, phi=F, mu=F)), check.attributes = FALSE)
  expect_equal(resC3_mathematica$PREDICT$PHI,   unname(predict(resC3, data.frame("depth"=0:150), alpha=F, phi=T, mu=F)), check.attributes = FALSE)
  expect_equal(resC3_mathematica$PREDICT$MU,    unname(predict(resC3, data.frame("depth"=0:150), alpha=F, phi=F, mu=T)), check.attributes = FALSE)

  # BEWARE! UGLY CODE AHEAD!
  conf_ints <- confint(resC3, level = c(.99, .95))
  conf_mat <- matrix(NA, 9L, 4L); rowz <- rep(1:3, 3L); listz <- rep(1:3, each=3L)
  for(i in 1:9){
    conf_mat[i,] <- c(conf_ints$ci[[2L]][[listz[i]]][rowz[i], 1L], conf_ints$ci[[1L]][[listz[i]]][rowz[i], ], conf_ints$ci[[2L]][[listz[i]]][rowz[i], 2L])
  }
  conf_ints <- unname(cbind(conf_mat[,1:2], unlist(conf_ints$coefficients), conf_mat[,3:4]))
  expect_equal(resC3_mathematica$CONFINT, conf_ints)
})

#
#
#
#
#

context("  AL: Common - Custom Model ( Y ~ depth | depth + I(depth^2) | 1 )\n   ")

resC4 <- DirichReg(Y ~ depth | depth + I(depth^2) | 1, AL)

load("testdata/resC4.RData")

test_that("Model Estimation", {
  expect_equal(resC4_mathematica$MLE, resC4$logLik)
  expect_equal(resC4_mathematica$DEV, -2.0*resC4$logLik)
  expect_equal(resC4_mathematica$COEFS, unname(resC4$coefficients), check.attributes = FALSE)
  expect_equal(resC4_mathematica$SE, unname(resC4$se), check.attributes = FALSE)
  expect_equal(resC4_mathematica$Z, unname(resC4$coefficients / resC4$se), check.attributes = FALSE)
  expect_equal(resC4_mathematica$P, 2*pnorm(-abs(unname(resC4$coefficients / resC4$se))), check.attributes = FALSE)
  expect_equal(resC4_mathematica$HESSIAN, unname(resC4$hessian), check.attributes = FALSE)
  expect_equal(resC4_mathematica$VCOV, unname(resC4$vcov), check.attributes = FALSE)
})

test_that("Methods", {
  expect_equal(resC4_mathematica$NOBS, nobs(resC4), check.attributes = FALSE)
  expect_equal(resC4_mathematica$MLE, unclass(logLik(resC4)), check.attributes = FALSE)
  expect_equal(resC4_mathematica$NPAR, attributes(logLik(resC4))$df, check.attributes = FALSE)
  expect_equal(resC4_mathematica$AIC, AIC(resC4), check.attributes = FALSE)
  expect_equal(resC4_mathematica$BIC, BIC(resC4), check.attributes = FALSE)
  expect_equal(resC4_mathematica$COEFS, unlist(coef(resC4)), check.attributes = FALSE)
  expect_equal(resC4_mathematica$VCOV, vcov(resC4), check.attributes = FALSE)

#  expect_equal(resC4_mathematica$PREDICT$ALPHA, unname(fitted(resC4, alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
#  expect_equal(resC4_mathematica$PREDICT$PHI,   unname(fitted(resC4, alpha=F, phi=T, mu=F))[1], check.attributes = FALSE)
#  expect_equal(resC4_mathematica$PREDICT$MU,    unname(fitted(resC4, alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  expect_equal(resC4_mathematica$PREDICT$ALPHA, unname(predict(resC4, data.frame("depth"=0:150), alpha=T, phi=F, mu=F)), check.attributes = FALSE)
  expect_equal(resC4_mathematica$PREDICT$PHI,   unname(predict(resC4, data.frame("depth"=0:150), alpha=F, phi=T, mu=F)), check.attributes = FALSE)
  expect_equal(resC4_mathematica$PREDICT$MU,    unname(predict(resC4, data.frame("depth"=0:150), alpha=F, phi=F, mu=T)), check.attributes = FALSE)

  # BEWARE! UGLY CODE AHEAD!
  conf_ints <- confint(resC4, level = c(.99, .95))
  conf_mat <- matrix(NA, 6L, 4L); rowz <- c(1,2,1,2,3,1); listz <- rep(1:3, c(2, 3, 1))
  for(i in 1:6){
    conf_mat[i,] <- c(conf_ints$ci[[2L]][[listz[i]]][rowz[i], 1L], conf_ints$ci[[1L]][[listz[i]]][rowz[i], ], conf_ints$ci[[2L]][[listz[i]]][rowz[i], 2L])
  }
  conf_ints <- unname(cbind(conf_mat[,1:2], unlist(conf_ints$coefficients), conf_mat[,3:4]))
  expect_equal(resC4_mathematica$CONFINT, conf_ints)
})
