# last modified 2013-12-19 by J. Fox


summary.objectiveML <- function(object, digits=getOption("digits"), conf.level=.90, robust=FALSE, analytic.se=object$t <= 500, 
                                fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"), ...) {
    fit.indices <- if (is.null(fit.indices)) ""
    else {
        if (missing(fit.indices)){
            if (is.null(opt <- getOption("fit.indices"))) c("AIC", "BIC") else opt
        }
        else match.arg(fit.indices, several.ok=TRUE)
    }
	vcov <- vcov(object, robust=robust, analytic=analytic.se)
	if (any(is.na(vcov))) stop("coefficient covariances cannot be computed")
	norm.res <- normalizedResiduals(object)
	se <- sqrt(diag(vcov))
	z <- object$coeff/se
	n.fix <- object$n.fix
	n <- object$n
	t <- object$t
	S <- object$S
	C <- object$C
	N <- object$N
	df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
	dfNull <- n*(n - 1)/2
	invC <- solve(C)
	CSC <- invC %*% (S - C)
	CSC <- CSC %*% CSC
	CS <- invC %*% S
	CS <- CS %*% CS
	chisqNull <- chisqNull(object)
	chisq <- object$criterion * (N - (!object$raw))
	GFI <- if (!"GFI" %in% fit.indices) NA else if (!object$raw) 1 - sum(diag(CSC))/sum(diag(CS)) else NA
    Rsq <- if (!object$raw) Rsq(object) else NA
	if ((!object$raw) && df > 0){
		AGFI <- if (!"AGFI" %in% fit.indices) NA else 1 - (n*(n + 1)/(2*df))*(1 - GFI)
		NFI <- if (!"NFI" %in% fit.indices) NA else (chisqNull - chisq)/chisqNull
		NNFI <- if (!"NNFI" %in% fit.indices) NA else (chisqNull/dfNull - chisq/df)/(chisqNull/dfNull - 1)
		L1 <- max(chisq - df, 0)
		L0 <- max(L1, chisqNull - dfNull)
		CFI <- if (!"CFI" %in% fit.indices) NA else 1 - L1/L0
        RNI <- if (!"RNI" %in% fit.indices) NA else 1 - (chisq - df)/(chisqNull - dfNull)
        IFI <- if (!"IFI" %in% fit.indices) NA else (chisqNull - chisq)/(chisqNull - df)
        if (!"RMSEA" %in% fit.indices) {
            RMSEA <- NA
            }
		else {
		    RMSEA <- sqrt(max(object$criterion/df - 1/(N - (!object$raw)), 0))
		    tail <- (1 - conf.level)/2 
		    max <- N
		    while (max > 1){
		        res <- optimize(function(lam) (tail - pchisq(chisq, df, ncp=lam))^2, interval=c(0, max))
		        if (is.na(res$objective) || res$objective < 0){
		            max <- 0
		            warning("cannot find upper bound of RMSEA")
		            break
		        }				
		        if (sqrt(res$objective) < tail/100) break
		        max <- max/2
		    }
		    lam.U <- if (max <= 1) NA else res$minimum
		    max <- max(max, 1)
		    while (max > 1){
		        res <- optimize(function(lam) (1 - tail - pchisq(chisq, df, ncp=lam))^2, interval=c(0, max))
		        if (sqrt(res$objective) < tail/100) break
		        max <- max/2
		        if (is.na(res$objective) || res$objective < 0){
		            max <- 0
		            warning("cannot find lower bound of RMSEA")
		            break
		        }				
		    }
		    lam.L <- if (max <= 1) NA else res$minimum
		    RMSEA.U <- sqrt(lam.U/((N - (!object$raw))*df))
		    RMSEA.L <- sqrt(lam.L/((N - (!object$raw))*df))
		}
	}
	else RMSEA.U <- RMSEA.L <- RMSEA <- NFI <- NNFI <- CFI <- AGFI <- RNI <- IFI <- NA
	if (!is.na(RMSEA)) RMSEA <- c(RMSEA, RMSEA.L, RMSEA.U, conf.level)
	if (!is.null(object$coeff)){
		var.names <- rownames(object$A)
		ram <- object$ram[object$par.posn, , drop=FALSE]
		par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
				var.names[ram[,3]])
		coeff <- data.frame(object$coeff, se, z, 2*pnorm(abs(z), lower.tail=FALSE), par.code)
		names(coeff) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)", " ")
		row.names(coeff) <- names(object$coeff)
	}
	else coeff <- NULL
	AIC <- if (!"AIC" %in% fit.indices) NULL else AIC(object)
	AICc <- if (!"AICc" %in% fit.indices) NULL else AICc(object)
	BIC <- if (!"BIC" %in% fit.indices) NULL else BIC(object)
	CAIC <- if (!"CAIC" %in% fit.indices) NULL else CAIC(object)
	SRMR <- if (!"SRMR" %in% fit.indices) NA else sqrt(sum(standardizedResiduals(object)^2 * 
							upper.tri(diag(n), diag=TRUE))/(n*(n + 1)/2))
	if (robust) { 
		chisq.adjusted <- object$adj.obj$chisq.scaled
		chisqNull.adjusted <- chisqNull/object$adj.obj$c 
		NFI.adjusted <- if (!"NFI" %in% fit.indices) NULL else (chisqNull.adjusted - chisq)/chisqNull.adjusted
		NNFI.adjusted <- if (!"NNFI" %in% fit.indices) NULL else (chisqNull.adjusted/dfNull - chisq.adjusted/df)/(chisqNull.adjusted/dfNull - 1)
		L1 <- max(chisq.adjusted - df, 0)
		L0 <- max(L1, chisqNull.adjusted - dfNull)
		CFI.adjusted <- if (!"CFI" %in% fit.indices) NULL else 1 - L1/L0
        RNI.adjusted <- if (!"RNI" %in% fit.indices) NULL else 1 - (chisq.adjusted - df)/(chisqNull.adjusted - dfNull)
        IFI.adjusted <- if (!"IFI" %in% fit.indices) NULL else (chisqNull.adjusted - chisq.adjusted)/(chisqNull.adjusted - df)
	}
	else{
		chisq.adjusted <- chisqNull.adjusted <- NFI.adjusted <- NNFI.adjusted <- CFI.adjusted <- RNI.adjusted <- IFI.adjusted <- NULL
	}
	ans <- list(chisq=chisq, df=df, chisqNull=chisqNull, dfNull=dfNull,
			GFI=GFI, AGFI=AGFI, RMSEA=RMSEA, NFI=NFI, NNFI=NNFI, CFI=CFI, RNI=RNI, IFI=IFI, BIC=BIC, SRMR=SRMR, 
			AIC=AIC, AICc=AICc, CAIC=CAIC, Rsq=Rsq,
			chisq.adjusted=chisq.adjusted, chisqNull.adjusted=chisqNull.adjusted, NFI.adjusted=NFI.adjusted,
			NNFI.adjusted=NNFI.adjusted, CFI.adjusted=CFI.adjusted, RNI.adjusted=RNI.adjusted, IFI.adjusted=IFI.adjusted,
			norm.res=norm.res, coeff=coeff, digits=digits, 
			iterations=object$iterations, aliased=object$aliased, raw=object$raw,
			robust=robust, robust.vcov=object$robust.vcov, adj.obj=object$adj.obj)
	class(ans) <- "summary.objectiveML"
	ans
}

print.summary.objectiveML <- function(x, digits=getOption("digits"), ...){
	old.digits <- options(digits=digits)
	on.exit(options(old.digits))
	if (x$raw) cat("\nModel fit to raw moment matrix.\n")	
	if (x$robust && !is.null(x$robust.vcov)){
		cat("\nSatorra-Bentler Corrected Fit Statistics:\n")
		cat("\n Corrected Model Chisquare = ", x$chisq.adjusted, "  Df = ", x$df, 
				"Pr(>Chisq) =", if (x$df > 0) pchisq(x$chisq.adjusted, x$df, lower.tail=FALSE)
						else NA)
		if (!x$raw) {		
			cat("\n Corrected Chisquare (null model) = ", x$chisqNull.adjusted,  "  Df = ", x$dfNull)
		}
		if (x$df > 0 && !x$raw){
			if (!is.null(x$NFI.adjusted)) cat("\n Corrected Bentler-Bonett NFI = ", x$NFI.adjusted)
			if (!is.null(x$NNFI.adjusted)) cat("\n Corrected Tucker-Lewis NNFI = ", x$NNFI.adjusted)
			if (!is.null(x$CFI.adjusted)) cat("\n Corrected Bentler CFI = ", x$CFI.adjusted)
            if (!is.null(x$RNI.adjusted)) cat("\n Corrected Bentler RNI = ", x$RNI.adjusted)
            if (!is.null(x$IFI.adjusted)) cat("\n Corrected Bollen IFI = ", x$IFI.adjusted)
		}
		cat("\n\nUncorrected Fit Statistics:\n")
		x$coeff[,2] <- sqrt(diag(x$robust.vcov))
		x$coeff[,3] <- x$coeff[,1]/x$coeff[,2]
		x$coeff[,4] <- 2*pnorm(abs(x$coeff[,3]), lower.tail=FALSE)
		colnames(x$coeff)[2] <- "Corrected SE"
	}
	if (!is.null(x$chisq)) cat("\n Model Chisquare = ", x$chisq, "  Df = ", x$df, 
				"Pr(>Chisq) =", if (x$df > 0) pchisq(x$chisq, x$df, lower.tail=FALSE)
						else NA)
	else if (!is.null(x$logLik)) cat("\n Model log-likelihood = ", x$logLik, "  Df = ", x$df, "\n")
	if (!x$raw) {		
		if ((!is.null(x$chisqNULL)) && (!is.na(x$chisqNULL))) cat("\n Chisquare (null model) = ", x$chisqNull,  "  Df = ", x$dfNull)
		if (!is.na(x$GFI)) cat("\n Goodness-of-fit index = ", x$GFI)
	}
	if (x$df > 0 && !x$raw){
		if (!is.na(x$AGFI)) cat("\n Adjusted goodness-of-fit index = ", x$AGFI)
		if (length(x$RMSEA) > 1 || !is.na(x$RMSEA)) cat("\n RMSEA index =  ", x$RMSEA[1],
					"   ", 100*x$RMSEA[4], "% CI: (", x$RMSEA[2], ", ", x$RMSEA[3],")", sep="")
		if (!is.na(x$NFI)) cat("\n Bentler-Bonett NFI = ", x$NFI)
		if (!is.na(x$NNFI)) cat("\n Tucker-Lewis NNFI = ", x$NNFI)
		if (!is.na(x$CFI)) cat("\n Bentler CFI = ", x$CFI)
        if (!is.na(x$RNI)) cat("\n Bentler RNI = ", x$RNI)
        if (!is.na(x$IFI)) cat("\n Bollen IFI = ", x$IFI)
		if (!is.na(x$SRMR)) cat("\n SRMR = ", x$SRMR)
	}
	if (!is.null(x$AIC) && !is.na(x$AIC)) cat("\n AIC = ", x$AIC)
	if (!is.null(x$AICc) && !is.na(x$AICc)) cat("\n AICc = ", x$AICc)
	if (!is.null(x$BIC) && !is.na(x$BIC)) cat("\n BIC = ", x$BIC)
	if (!is.null(x$CAIC) && !is.na(x$CAIC)) cat("\n CAIC = ", x$CAIC)
	if (length(x$norm.res) > 1 || !is.na(x$norm.res)){
		cat("\n\n Normalized Residuals\n")
		print(summary(as.vector(x$norm.res)))
	}
	if (!is.na(x$Rsq[1])){
		cat("\n R-square for Endogenous Variables\n")
		print(round(x$Rsq, 4))
	}
	if (!is.null(x$coeff)){
		if (x$robust && !is.null(x$robust.vcov)) cat("\n Parameter Estimates (with Robust Standard Errors)\n") else cat("\n Parameter Estimates\n")
		print(x$coeff, right=FALSE, digits=digits)
		if (!is.na(x$iterations)) cat("\n Iterations = ", x$iterations, "\n")
		if (!is.null(x$aliased)) cat("\n Aliased parameters:", x$aliased, "\n")
	}
	invisible(x)
}

summary.objectiveGLS <- function(object, digits=getOption("digits"), conf.level=.90, 
                                 fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR"), ...){
    fit.indices <- if (missing(fit.indices)){
        getOption("fit.indices")
    }
    else if (fit.indices[1] != "") match.arg(fit.indices, several.ok=TRUE)
	summary <-  summary.objectiveML(object, digits=digits, conf.level=conf.level, analytic.se=FALSE, fit.indices=fit.indices, ...)
	S <- object$S
	Sinv <- solve(S)
	C <- object$C
	SinvSmC <- Sinv %*% (S - C)
	SinvS <- Sinv %*% S
	n <- object$n
	if ("GFI" %in% fit.indices) summary$GFI <- 1 - sum(diag(SinvSmC %*% SinvSmC))/sum(diag(SinvS %*% SinvS))
	if ("AGFI" %in% fit.indices) summary$AGFI <-  1 - (n*(n + 1)/(2*summary$df))*(1 - summary$GFI)
	summary
}

deviance.objectiveML <- function(object, ...) object$criterion * (object$N - (!object$raw))

df.residual.sem <- function(object, ...) {
	n.fix <- object$n.fix
	n <- object$n
	t <- object$t
	n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
}

Rsq <- function(model){
	A <- model$A
	P <- model$P
	IAinv <- solve(diag(nrow(A)) - A)
	C <- IAinv %*% P %*% t(IAinv)
	R2 <- 1 - diag(P)/diag(C)
	R2 <- R2[classifyVariables(model$semmod)$endogenous]
	R2
}
