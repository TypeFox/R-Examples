mfp.fit <- function(x, y, cox, gauss, df, scaling, alpha, select, verbose = TRUE, xnames = NULL, maxits = 20, ...)
{
    int <- as.numeric(!cox) # intercept
    nx <- ncol(x) - int
    nobs <- nrow(x)
    x.names <- dimnames(x)[[2]][int + seq(nx)]
#    if(sum((df == 1 | df == 2 | df == 4), na.rm=TRUE) != nx)
    if(sum(df %in% c(1,2,4), na.rm=TRUE) != nx)
        stop("df is invalid")
    if(sum((alpha > 0 & alpha <= 1), na.rm=TRUE) != nx)
        stop("alpha is invalid")
    if(sum((select > 0 & select <= 1), na.rm=TRUE) != nx) stop(
            "select is invalid")    
	if(is.null(xnames)) xnames <- x.names
#
# Step 1: Order variables by LR test
#  via one-step backward selection
#
    x.order <- fp.order(x, y, cox, gauss, xnames, ...)
    x <- x[, c(int, int + x.order$order), drop=FALSE]
    x.names <- x.names[x.order$order]
    df <- df[x.order$order]
    alpha <- alpha[x.order$order]
    select <- select[x.order$order]
    scaling <- scaling[x.order$order]
#
# Set up powers & working matrix
#
    pwrs.mx <- matrix(c(rep(1, nx), rep(NA, nx)), ncol = 2, byrow=FALSE, 
        dimnames = list(x.names, c("power1", "power2")))    
    # NA = no power trafo
    pvals.mx <- matrix(rep(NA, nx * 6), ncol = 6, dimnames = list(x.names,
        c("p.null", "p.lin", "p.FP", "power2", "power4.1", "power4.2")))
    scale.mx <- matrix(0:1, ncol = 2, nrow = nx, byrow = TRUE, dimnames = list(x.names, c(
        "shift", "scale")))
    df.work <- rep(1, nx)
    if(nx > 1)
        pwrs.comp <- matrix(ncol = 2, nrow = nx * (nx - 1))
#
# x.work: working matrix with doubled no of columns for input vars.
# initially second column entries are set to 0 
#
    x.work <- matrix(0, nrow = nobs, ncol = int + 2 * nx)
    if(!cox)
        x.work[, 1] <- x[, 1]
    x.work[, 2 * seq(nx) - 1 + int] <- x[, seq(nx) + int]
    x.work.names <- c(paste(rep(x.names, each = 2), ".", rep(1:2, nx), sep
         = ""))
    if(cox)
        dimnames(x.work) <- list(1:nobs, x.work.names)
    else dimnames(x.work) <- list(1:nobs, c("Intercept", x.work.names))
#
    pwrs.stable <- FALSE    
#
# Step 2: Backfitting loop
#
    its <- 0
	if(verbose) {
		pos <- c(1, 3, 5) # output formatting
		fp.out("Variable", "Deviance", "Power(s)", pos = pos)
		cat("\n", rep("-", 48), sep = "")
	}
    while(!pwrs.stable & its < maxits) {
        its <- its + 1
        j <- 0
#        while(j < nx & !pwrs.stable) {
        while(j < nx) {
            j <- j + 1
# Check convergence 
            if(nx > 1) {
                num <- (j - 1) * (nx - 1) + seq(nx - 1)
                if(its > 1) {
					mx.old <- pwrs.comp[num,, drop=FALSE]
					mx.new <- pwrs.mx[-j,, drop=FALSE]
#
					pwrs.stable <- mx.com(mx.old, mx.new)
                }
                pwrs.comp[num,] <- pwrs.mx[-j,, drop=FALSE]
                dfr <- nobs - sum(df.work[-j]) - int
            }
            else dfr <- nobs - int
            if(!pwrs.stable) {
#
# Step 2a: Set up matrices and fit and find best single FP. Start with best drop1 variable
#
                num <- 2 * (j - 1) + int + seq(2)    # possible new positions of power-trafo var
                xj <- x[, j + int]
                if(its == 1) {
					if(df[j]>1) {
						if(min(xj) <= 0 & scaling[j] == FALSE) 
							stop(paste("\nSome data values in '",x.names[[j]],"' are less or equal zero. To use fp all data must be > 0.\nApply scale=TRUE or use appropriately shifted data instead.",sep=""))
						xj.transform <- fp.scale(xj, scaling[j])
						scale.mx[j, 1] <- xj.transform$shift
						scale.mx[j, 2] <- xj.transform$scale
					}
                }
                fitj <- fp.fit(cbind(xj, x.work[,  -num, drop=FALSE]), y, df[j], dfr, cox, gauss, scale.mx[j, 1], scale.mx[j, 2], ...)
                res <- fp.sel(fitj, alpha[j], select[j])
                best.fitj <- res$results
                fit.fitj <- res$fit
                pwrsj <- best.fitj$pwrs
                pwrs.mx[j,  ] <- pwrsj
                df.work[j] <- best.fitj$df
                x.work[, num] <- fp.gen(xj, pwrsj, scale.mx[j, 1], scale.mx[j, 2])   
# Verbose output
				if(verbose) {
					if(j==1) cat("\nCycle", its)
					pos <- c(1, 3, 5)   
					namej <- x.names[[j]]
					fp.out(namej," ", " ", pos = pos)
					fp.out("        ", round(fit.fitj$dev0, 3), " ", pos = pos)
					fp.out("        ", round(fit.fitj$dev1, 3), "1", pos = pos)
					fp.out("        ", round(fit.fitj$dev2, 3), fit.fitj$pwr2, pos = pos)
					fp.out("        ", round(fit.fitj$dev4, 3), fit.fitj$pwr4, pos = pos)
					cat("\n")
				}
#
				if(df[j] < 4)
					pvals.mx[j,  ] <- c(best.fitj$p.null, best.fitj$p.lin, best.fitj$p.FP, fitj$pwr2, NA, NA)
				else pvals.mx[j,  ] <- c(best.fitj$p.null, best.fitj$p.lin, best.fitj$p.FP, fitj$pwr2, fitj$pwr4)
			} # end of "if(!powers.stable)"
		}   # end of "while(j < nx)"
		if(nx == 1)
            pwrs.stable <- 1
    } # end of "while(!pwrs.stable)"
    if(verbose) cat("\n")   # end
	if(its==maxits) stop(paste("No convergence within",maxits,"iterations."))
#
# Set up final matrix
#
    num <- int + seq(2 * nx)[is.na(as.vector(t(pwrs.mx)))]
    if(sum(num, na.rm=TRUE) > 0)
        x.work <- x.work[,  - num, drop=FALSE]
    if(cox) {
        if(exists("coxph.fit")) fitter <- get("coxph.fit")
        else fitter <- getFromNamespace("coxph.fit","survival")
    }
    else {
        fitter <- get("glm.fit")
    }
    fit <- fitter(x.work, y, ...)
    fit$x <- x.work
    fit$powers <- pwrs.mx
    fit$pvalues <- pvals.mx
    fit$scale <- scale.mx
    fit$df.initial <- matrix(df, dimnames = list(x.names, "df.initial"))
    fit$df.final <- matrix(df.work, dimnames = list(x.names, "df.final"))
    fit$dev <- best.fitj$dev
    fit$dev.lin <- x.order$dev[2]
    fit$df.lin <- x.order$df[2]
    fit$dev.null <- x.order$dev[1]
    fit$df.null <- x.order$df[1]
	if(cox) fit$method <- "efron"
#
# Final output
#
    power1 <- pwrs.mx[, 1]
    power2 <- pwrs.mx[, 2]
    power1[is.na(power1)] <- "."
    power2[is.na(power2)] <- "."
    fit$fptable <- data.frame(df.initial = df, select, alpha, df.final = df.work,
        power1, power2, row.names = x.names)
if(verbose) {
    cat("\nTansformation\n")
    print(fit$scale)
    cat("\nFractional polynomials\n")
    print(fit$fptable)
    cat("\n")
}
fit
}
