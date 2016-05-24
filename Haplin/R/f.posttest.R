f.posttest <- function(coef_, cov_, test){
##
## EXTRACT COEFFICIENTS TO BE TESTED
## BUILD CONTRAST MATRIX FOR INTERACTIONS
## DO THE CHI-SQUARED TEST
##
.vis <- F
#
.coef <- coef_
.cov <- cov_
#
## SELECT COEFFICIENTS TO BE TESTED
f.vis(.coef <- lapply(.coef, function(x)x[test, , drop = F]), vis = .vis)
f.vis(.cov <- lapply(.cov, function(x) x[test, test, drop = F]), vis = .vis)
#
## RESHAPE COEFFICIENTS AND COVARIANCE MATR. INTO FULL SIZE
.n.pars <- length(.coef[[1]])
.l <- length(.coef)
f.vis(.coef.vec <- unlist(.coef), vis = .vis)
f.vis(.cov.mat <- f.bdiag(.cov), vis = .vis)
#
## BUILD CONTRAST MATRIX
.A <- f.post.contrasts(test.type = "interaction", n.res = .l, n.pars = .n.pars)
#
## DO CHI-SQUARED TEST
.chisq.res <- f.post.chisq(coeff = .coef.vec, covar = .cov.mat, contrast.mat = .A)
#
##
return(.chisq.res)
}
