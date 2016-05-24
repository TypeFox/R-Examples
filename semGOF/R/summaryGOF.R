chisqNull <- function(object){
	UseMethod("chisqNull")
}

chisqNull.objectiveML <- function(object){
	chisq <- if (!object$raw) {
			S <- object$S
			CC <- diag(diag(S))
			(object$N - 1) * 
				(sum(diag(S %*% solve(CC))) + log(det(CC)) - log(det(S)) - object$n)
		}
		else NULL
	chisq
}  

summaryGOF <- function(object, digits=5, ...) {
	vcov <- vcov(object)
	if (any(is.na(vcov))) stop("coefficient covariances cannot be computed")
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
	
        if (!object$raw){
        ICOMP <- - log(chisq) + 0.5 * (t + n.fix*(n.fix + 1)/2) * log((sum(diag(C)))/(t + n.fix*(n.fix + 1)/2)) - 0.5 * log(det(C))
        Fml <- log(det(C)) - log(det(S)) + sum(diag(invC %*% S)) - n
        }

	if ((!object$raw) && df > 0){
	RNI <- 1 - (chisq - df)/(chisqNull - dfNull)
        IFI <- (chisqNull - chisq)/(chisqNull - df)
        chisq.df <- ((N - 1)/df) * (chisq/(N - 1))
        CN <- ((1.96 + sqrt(2 * df - 1))^2/(2 * chisq/(N-1))) + 1
        Gamma.hat <-  n/(n + 2 * ((chisq - df)/(N - 1))) 
        BL86 <- ((chisqNull/dfNull) - (chisq/df))/(chisqNull/dfNull)
        W <- chisq/df
		
	}
	lambda <- chisq - df
        d <- lambda/N
        Mc <- exp(-0.5 * d)
        CAK <- (chisq/(N - 1)) + 2 * n/N
        CSK <- (chisq/(N - 1)) + (n * log(N))/N
        ECVI <- chisq/(N - 1) + 2 * ((t + n.fix*(n.fix + 1)/2)/(N - 1))

	ob <- list(ICOMP=ICOMP, Fml=Fml, RNI=RNI, IFI=IFI, chisq.df=chisq.df, CN=CN, Gamma.hat=Gamma.hat, BL86=BL86, W=W, d=d, Mc=Mc, CAK=CAK, CSK=CSK, ECVI=ECVI, digits=digits, aliased=object$aliased, raw=object$raw, df=df)
	class(ob) <- "summaryGOF"
	ob
}

print.summaryGOF <- function(x, ...){
	old.digits <- options(digits=x$digits)
	on.exit(options(old.digits))
        cat("\n Goodness-of-Fit indexes of structural equation models for 'sem' package\n")
	if (x$raw) cat("\nModel fit to raw moment matrix.\n")	
	
	if (!x$raw) {		
		cat("\n ICOMP = ", x$ICOMP)
                cat("\n Fml = ", x$Fml)
	}
	if (x$df > 0 && !x$raw){
		cat("\n RNI = ", x$RNI)
                cat("\n IFI = ", x$IFI)
                cat("\n chisq.df = ", x$chisq.df)
                cat("\n CN = ", x$CN)
                cat("\n Gamma.hat = ", x$Gamma.hat)
                cat("\n BL86 = ", x$BL86)
                cat("\n W = ", x$W)
	
	}
      cat("\n d = ", x$d)
      cat("\n Mc = ", x$Mc)
      cat("\n CAK = ", x$CAK)
      cat("\n CSK = ", x$CSK)
      cat("\n ECVI = ", x$ECVI, "\n")
	
	invisible(x)
}


