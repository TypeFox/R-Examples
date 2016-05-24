#TO DO:
#anova(resC1)
#residuals(resC1)
#update(resC1)

tol3 <- .Machine$double.eps^(1/3)

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

cat("\n  --- Checks: Alternative Model ------------------------------------------------")

##                                                                             #
##                                                                            ##
##                                                                             #
##                                                                             #
##                                                                             #

context("  AL: Alternative - Null Model ( Y ~ 1 ), base = 1\n   ")

resA1_1 <- DirichReg(Y ~ 1, AL, model = "alternative", base = 1L)

load("testdata/resA1_1.RData")

test_that("Model Estimation", {
  expect_equal(resA1_1_mathematica$MLE, resA1_1$logLik)
  expect_equal(resA1_1_mathematica$DEV, -2.0*resA1_1$logLik)
  expect_equal(resA1_1_mathematica$COEFS, unname(resA1_1$coefficients), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$SE, unname(resA1_1$se), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$Z, unname(resA1_1$coefficients / resA1_1$se), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$P, 2*pnorm(-abs(unname(resA1_1$coefficients / resA1_1$se))), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$HESSIAN, unname(resA1_1$hessian), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$VCOV, unname(resA1_1$vcov), check.attributes = FALSE)
})

test_that("Methods", {
  expect_equal(resA1_1_mathematica$NOBS, nobs(resA1_1), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$MLE, unclass(logLik(resA1_1)), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$NPAR, attributes(logLik(resA1_1))$df, check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$AIC, AIC(resA1_1), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$BIC, BIC(resA1_1), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$COEFS, unlist(coef(resA1_1)), check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$VCOV, vcov(resA1_1), check.attributes = FALSE)

  expect_equal(resA1_1_mathematica$PREDICT$ALPHA, unname(fitted(resA1_1, alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$PREDICT$PHI,   unname(fitted(resA1_1, alpha=F, phi=T, mu=F))[1], check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$PREDICT$MU,    unname(fitted(resA1_1, alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  expect_equal(resA1_1_mathematica$PREDICT$ALPHA, unname(predict(resA1_1, data.frame("depth"=0), alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$PREDICT$PHI,   unname(predict(resA1_1, data.frame("depth"=0), alpha=F, phi=T, mu=F))[1,], check.attributes = FALSE)
  expect_equal(resA1_1_mathematica$PREDICT$MU,    unname(predict(resA1_1, data.frame("depth"=0), alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  conf_ints <- confint(resA1_1, level = c(.99, .95))
  conf_ints$coefficients[[1L]][4L] <- conf_ints$coefficients[[2L]]
  conf_ints <- lapply(c(2L, 3L, 4L), function(listelement){ sort(unlist(lapply(c(conf_ints$coefficients[1L], conf_ints$ci), `[[`, listelement))) })
  conf_ints <- unname(t(matrix(unlist(conf_ints), 5)))
  expect_equal(resA1_1_mathematica$CONFINT, conf_ints)
})

##                                                                          ###
##                                                                         #   #
##                                                                            #
##                                                                          # 
##                                                                         #####

context("  AL: Alternative - Linear/Constant Model ( Y ~ depth | 1 ), base = 1\n   ")

resA2_1 <- DirichReg(Y ~ depth | 1, AL, model = "alternative", base = 1L)

load("testdata/resA2_1.RData")

test_that("Model Estimation", {
  expect_equal(resA2_1_mathematica$MLE, resA2_1$logLik)
  expect_equal(resA2_1_mathematica$DEV, -2.0*resA2_1$logLik)
  expect_equal(resA2_1_mathematica$COEFS, unname(resA2_1$coefficients), check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$SE, unname(resA2_1$se), check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$Z, unname(resA2_1$coefficients / resA2_1$se), check.attributes = FALSE, tolerance = tol3)
  expect_equal(resA2_1_mathematica$P, 2*pnorm(-abs(unname(resA2_1$coefficients / resA2_1$se))), check.attributes = FALSE, tolerance = tol3)
  expect_equal(resA2_1_mathematica$HESSIAN, unname(resA2_1$hessian), check.attributes = FALSE, tolerance = tol3)
  expect_equal(resA2_1_mathematica$VCOV, unname(resA2_1$vcov), check.attributes = FALSE, tolerance = tol3)
})

test_that("Methods", {
  expect_equal(resA2_1_mathematica$NOBS, nobs(resA2_1), check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$MLE, unclass(logLik(resA2_1)), check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$NPAR, attributes(logLik(resA2_1))$df, check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$AIC, AIC(resA2_1), check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$BIC, BIC(resA2_1), check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$COEFS, unlist(coef(resA2_1)), check.attributes = FALSE)
  expect_equal(resA2_1_mathematica$VCOV, vcov(resA2_1), check.attributes = FALSE, tolerance = tol3)

#  expect_equal(resA2_1_mathematica$PREDICT$ALPHA, unname(fitted(resA2_1, alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
#  expect_equal(resA2_1_mathematica$PREDICT$PHI,   unname(fitted(resA2_1, alpha=F, phi=T, mu=F))[1], check.attributes = FALSE)
#  expect_equal(resA2_1_mathematica$PREDICT$MU,    unname(fitted(resA2_1, alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  expect_equal(resA2_1_mathematica$PREDICT$ALPHA, unname(predict(resA2_1, data.frame("depth"=0:150), alpha=T, phi=F, mu=F)), check.attributes = FALSE, tolerance = tol3)
  expect_equal(resA2_1_mathematica$PREDICT$PHI[1L],   unname(predict(resA2_1, data.frame("depth"=0:150), alpha=F, phi=T, mu=F))[1L], check.attributes = FALSE, tolerance = tol3)
  expect_equal(resA2_1_mathematica$PREDICT$MU,    unname(predict(resA2_1, data.frame("depth"=0:150), alpha=F, phi=F, mu=T)), check.attributes = FALSE)

  conf_ints <- confint(resA2_1, level = c(.99, .95))
  conf_ints$coefficients[[1L]][4L] <- conf_ints$coefficients[[2L]]
  conf_ints <- lapply(c(2L, 3L, 4L), function(listelement){ sort(unlist(lapply(c(conf_ints$coefficients[1L], conf_ints$ci), `[[`, listelement))) })
  conf_ints <- unname(t(matrix(unlist(conf_ints), 5)))
  expect_equal(resA2_1_mathematica$CONFINT, conf_ints)
})

##                                                                         ####
##                                                                             #
##                                                                           ##
##                                                                             #
##                                                                         ####

context("  AL: Alternative - Linear/Constant Model ( Y ~ 1 | depth ), base = 2\n   ")

resA3_2 <- DirichReg(Y ~ 1 | depth, AL, model = "alternative", base = 2L)

resA3_2_mathematica <- list(
  MLE     =  51.71895090427818151113,
  DEV     = -103.4379018085563630223,
  COEFS   = c(-1.92806496209330522141400, -0.13956410084104567041840, -0.07823484089230527185368, 0.04951503424822068330203),
  SE      = c(0.1762725713617107560414, 0.07280101211689304403532, 0.2818904517845701991150, 0.006980898589083690422078),
  Z       = c(-10.93797490556214811640, -1.917062644911507304984, -0.2775363280204143549263, 7.092931320568002149793),
  P       = c(7.587578202342218574969e-28, 0.05522997357879426516548, 0.7813683152144271216545, 1.313006280096710930199e-12),
  GRAD    = c(8.603474621174026685438e-31, 1.275594390293843165233e-30, 2.552444010017015860838e-31, 7.883692762433664995613e-28),
  HESSIAN = matrix(c(-89.23891691615910553897, 41.63244343038315899138, -22.52830663115145384000, -2658.038920499177560877, 41.63244343038315899138, -255.1197075452105633482, 8.794723038895463272885, 2278.390921785228498514, -22.52830663115145384000, 8.794723038895463272885, -51.03870689707925201174, -2219.570401094650987676, -2658.038920499177560877, 2278.390921785228498514, -2219.570401094650987676, -157673.8539428603750109), 4L),
  VCOV    = matrix(c(0.03107201941446941098236, -0.002278463636880385500344, 0.02605171520347625214022, -0.0009234599150738164961845, -0.002278463636880385500344, 0.005299987365244007816728, -0.007946771397992568261202, 0.0002268613711727447735950, 0.02605171520347625214022, -0.007946771397992568261202, 0.07946222680730909615495, -0.001672593932342348746343, -0.0009234599150738164961845, 0.0002268613711727447735950, -0.001672593932342348746343, 0.00004873294511107065961980), 4L),
  NOBS    = 39L,
  NPAR    = 4L,
  AIC     = -95.43790180855636302226,
  BIC     = -88.78365522403777731247,
  CONFINT = t(matrix(c(-2.382113016818714547595, -2.273552853424524888651, -1.928064962093305221414, -1.582577070762085554176, -1.474016907367895895232, -0.3270870811797573655124, -0.2822514626282201296519, -0.1395641008410456704184, 0.003123260946128788815203, 0.04795887949766602467567, -0.8043365269896397174919, -0.6307299739757875214791, -0.07823484089230527185368, 0.4742602921911769777717, 0.6478668452050291737845, 0.03153343109735573704559, 0.03583272443389017066214, 0.04951503424822068330203, 0.06319734406255119594191, 0.06749663739908562955846), 5L))#,
#  PREDICT = list(
#    ALPHA = c(1.021200212972279022750,2.318380244412313173929,1.298665556155223040582),
#    PHI   = 4.638246013539815237261,
#    MU    = c(0.2201694800127515751630,0.4998398613709953725498,0.2799906586162530522872)
#  )
)

test_that("Model Estimation", {
  expect_equal(resA3_2_mathematica$MLE, resA3_2$logLik)
  expect_equal(resA3_2_mathematica$DEV, -2.0*resA3_2$logLik)
  expect_equal(resA3_2_mathematica$COEFS, unname(resA3_2$coefficients), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$SE, unname(resA3_2$se), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$Z, unname(resA3_2$coefficients / resA3_2$se), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$P, 2*pnorm(-abs(unname(resA3_2$coefficients / resA3_2$se))), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$HESSIAN, unname(resA3_2$hessian), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$VCOV, unname(resA3_2$vcov), check.attributes = FALSE)
})

test_that("Methods", {
  expect_equal(resA3_2_mathematica$NOBS, nobs(resA3_2), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$MLE, unclass(logLik(resA3_2)), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$NPAR, attributes(logLik(resA3_2))$df, check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$AIC, AIC(resA3_2), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$BIC, BIC(resA3_2), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$COEFS, unlist(coef(resA3_2)), check.attributes = FALSE)
  expect_equal(resA3_2_mathematica$VCOV, vcov(resA3_2), check.attributes = FALSE)

#  expect_equal(resA3_2_mathematica$PREDICT$ALPHA, unname(fitted(resA3_2, alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
#  expect_equal(resA3_2_mathematica$PREDICT$PHI,   unname(fitted(resA3_2, alpha=F, phi=T, mu=F))[1], check.attributes = FALSE)
#  expect_equal(resA3_2_mathematica$PREDICT$MU,    unname(fitted(resA3_2, alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)
#
#  expect_equal(resA3_2_mathematica$PREDICT$ALPHA, unname(predict(resA3_2, data.frame("depth"=0), alpha=T, phi=F, mu=F))[1,], check.attributes = FALSE)
#  expect_equal(resA3_2_mathematica$PREDICT$PHI,   unname(predict(resA3_2, data.frame("depth"=0), alpha=F, phi=T, mu=F))[1,], check.attributes = FALSE)
#  expect_equal(resA3_2_mathematica$PREDICT$MU,    unname(predict(resA3_2, data.frame("depth"=0), alpha=F, phi=F, mu=T))[1,], check.attributes = FALSE)

  conf_ints <- confint(resA3_2, level = c(.99, .95))
  conf_ints <- unname(rbind(
    sort(c(conf_ints$coefficients$beta[[1L]], unlist(lapply(conf_ints$ci, `[`, 1L)))),
    sort(c(conf_ints$coefficients$beta[[3L]], unlist(lapply(conf_ints$ci, `[`, 3L)))),
    sort(c(conf_ints$coefficients$gamma[[1L]][[1L]], unlist(lapply(conf_ints$ci, function(lelment){ lelment$gamma[1L,] })))),
    sort(c(conf_ints$coefficients$gamma[[1L]][[2L]], unlist(lapply(conf_ints$ci, function(lelment){ lelment$gamma[2L,] }))))
  ))
  expect_equal(resA3_2_mathematica$CONFINT, conf_ints)
})
