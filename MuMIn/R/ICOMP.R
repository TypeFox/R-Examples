#ICOMP information criterion (I for informational and COMP for complexity)
#developed by Bozdogan (1988, 1990, 1993, 1994)

#Bozdogan, H (1990) On the information-based measure of covariance complexity
#and its application to the evaluation of multivariate linear models. \emph{Comm.
#Stat. Theory and Methods} 19:221-278
#Bozdogan, H. (2000). Akaike's Information Criteria and recent developments in
#Information Complexity. \emph{J. Math. Psych.} 44:62-91
#Bozdogan, H. and Haughton, D.M.A.  (1998) Information complexity criteria for
#regression models. \emph{Comp. Stat. & Data Analysis} 28:51-76


`ICOMP` <-
# function (object, ..., type = c("vcov", "r", "cv"), REML = NULL) {
function (object, ..., REML = NULL) {
    # type <- match.arg(type)
    type <- "vcov"
    ret <- sapply(list(object, ...), function(x) {
        ll <- if (!is.null(REML) && inherits(x, c("mer", "lme",
            "gls", "lm"))) logLik(x, REML = REML) else logLik(x)
        covmat <- vcov(x)
        k <- nrow(covmat) # attr(ll, "df")
        switch(type, vcov = {
            mat <- covmat
        }, r = {
            cov <- diag(diag(1 / covmat), nrow = nrow(covmat), ncol = ncol(covmat))
            mat <- sqrt(cov) %*% covmat %*% sqrt(cov)
        }, cv = {
            coefs <- coef(x)
            ncoef <- length(coefs)
            coefmat <- diag(1 / coefs, nrow = ncoef, ncol = ncoef)
            mat <- coefmat %*% covmat %*% coefmat
        })
        as.vector(-2 * c(ll) + k * log(sum(diag(mat)) / k) - log(det(mat)))
		# ICOMP = -2 * LL + k * log(tr(IFIM) / k) - log(det(IFIM))
		# ICOMP = -2 * LL + 2 * C * (sig(model))
		# where C is a complexity measure and sig(model) is the variance-
		# covariance matrix of the parameters estimated under the model.
    })
    if (length(ret) > 1L) {
        Call <- match.call()
        Call$type <- Call$REML <- NULL
        ret <- data.frame(ICOMP = ret, row.names = as.character(Call[-1L]))
    }
    return(ret)
}

# Bozdogans's CAICF (C denoting "consistent" and F denoting the use of the Fisher
# information matrix),

`CAICF` <-
function (object, ..., REML = NULL) {
    ret <- sapply(list(object, ...), function(x) {
        ll <- if (!is.null(REML) && inherits(x, c("mer", "lme",
            "gls", "lm"))) logLik(x, REML = REML) else logLik(x)
        covmat <- vcov(x)
        k <- attr(ll, "df")
		n <- attr(ll, "nobs")
		#log(det(solve(covmat))) == -log(det(covmat))
		-(2 * c(ll)) + (k * (log(n) + 2)) - log(det(covmat))
    })
	# I is the natural logarithm of the determinant of the estimated Fisher information matrix.
    if (length(ret) > 1L) {
        Call <- match.call()
        Call$type <- Call$REML <- NULL
        ret <- data.frame(CAICF = ret, row.names = make.unique(as.character(Call[-1L])))
    }
    return(ret)
}
