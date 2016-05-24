fp.fit <- function(X, Y, df, dfr, cox, gauss, shift, scale, ...)
{
#
# df defines DegreesOfFreedom of test and also degree of FP to test (df=1: linear, df=2: FP1, df=4: FP2)
#
    X <- X[,apply(X, 2, function(x) !all(x==0)), drop=FALSE]    # to avoid numerical problems
    x <- X[, 1]                         # x is first column of X	# induced by fp.gen if df=0
    Xcov <- X[, -1, drop=FALSE]         # all other covariates
    ncov <- ncol(Xcov)
    nobs <- nrow(X)
    pwrs <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
    npwrs <- length(pwrs)
    dev2 <- dev4 <- Inf                 # Deviance for FP1 and FP2 model
#
# Set up fitters
#
    if(cox) {
        if(exists("coxph.fit")) fitter <- get("coxph.fit")
        else fitter <- getFromNamespace("coxph.fit","survival")
        dv <- "loglik"
    }
    else {
        fitter <- get("glm.fit")
        dv <- "deviance"
    }
    n <- 1 + cox
#
# Null and linear models
#
	fit <- fitter(X, Y, ...)

	dispersion <- if (gauss) {
			if (fit$df.residual > 0) sum(fit$residuals^2)/fit$df.residual
			else Inf
			}
			else 1
#
    dev <- (-2)^cox * fit[[dv]]        # dev1: Full model resid deviance
    dev1 <- dev[n]
    if(!ncov)
        dev0 <- dev[1]
    else dev0 <- (-2)^cox * fitter(Xcov, Y, ...)[[dv]][n]  # dev0: Drop1 residual deviance
#
    if(df > 1) {
        for(i in 1:npwrs) {
#
# Find best single power transformation
#
            x.fp <- fp.gen(x, pwrs[i], shift, scale)
            dev <- (-2)^cox * fitter(cbind(x.fp, Xcov), Y, ...)[[dv]][n]
            if(!is.null(dev) & dev < dev2) {
                dev2 <- dev
                pwr2 <- pwrs[i]
            }
            if(df == 4) {
#
# Find best two power transformation
#
                x.fp <- fp.gen(x, c(pwrs[i], pwrs[i]), shift, scale)
                dev <- (-2)^cox * fitter(cbind(x.fp, Xcov), Y, ...)[[dv]][n]
                if(!is.null(dev) & dev < dev4) {
                  dev4 <- dev
                  pwr4 <- c(pwrs[i], pwrs[i])
                }
                j <- i + 1
                while(j <= npwrs) {
                  x.fp <- fp.gen(x, c(pwrs[i], pwrs[j]), shift, scale)
                  dev <- (-2)^cox * fitter(cbind(x.fp, Xcov), Y, ...)[[dv]][n]
                  if(!is.null(dev) & dev < dev4) {
                    dev4 <- dev
                    pwr4 <- c(pwrs[i], pwrs[j])
                  }
                  j <- j + 1
                }
            }
        }
    }
#
# Output
#   
    dev <- c(dev0, dev1, dev2, dev4)
    if(df < 4)
        dev[4] <- pwr4 <- NA
    if(df < 2)
        dev[3] <- pwr2 <- NA
#

# Compute -log Likelihood
 # if(gauss) {
 #	wt <- rep(1, nrow(X)) # actually use equal weights only
 #    dev <- sum(wt) * (log(dev/sum(wt) * 2 * pi) + 1) + 2
 # }
#
    fit <- list(pwr4 = pwr4, pwr2 = pwr2, dev4 = dev[4], dev2 = dev[3], 
        dev1 = dev[2], dev0 = dev[1], dispersion=dispersion, nobs = nobs, dfr = dfr, df = df, 
        gauss = gauss)
    return(fit)
}











