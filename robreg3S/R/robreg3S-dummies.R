.robreg3S_withDummies <- function(y, x, dummies, filter=TRUE, alpha=0.20, K=20, xi=0.01, ...){
  x <- as.matrix(x)
  dummies <- as.matrix(dummies)
  p <- ncol(x)
  pd <- ncol(dummies)
  n <- nrow(x)
  if( is.null(colnames(x)) ) 
    colnames(x) <- paste0("X",1:p)
  if( is.null(colnames(dummies)) ) 
    colnames(dummies) <- paste0("D",1:pd)

  ## filter step
  xf <- x
  if( filter )
    for(j in 1:p) xf[,j] <- .cfilter.iter(x[,j], alpha=alpha, miter=3)
  if( mean(rowSums(is.na(xf)) > 0) < xi ) xf <- x
  ## SM procedure
	bdummies <- suppressWarnings(as.numeric(coef(rlm( y ~ dummies, method="M"))))
	Bdummies <- apply(x, 2, function(x){
	  suppressWarnings(as.numeric(coef(rlm( x ~ dummies, method="M")	)))
	})
	yt <- y - cbind(1, dummies) %*% bdummies
	xt <- x - cbind(1, dummies) %*% Bdummies
	fit0 <- .robreg3S_noDummies(yt, xt, filter=filter, alpha=alpha, ...)
	b0 <- fit0$coef
	bdummies0 <- suppressWarnings(as.numeric(coef(rlm( y - cbind(1, fit0$ximp) %*% b0 ~ dummies, method="M")	)))
  for(i in 1:K){
	  yt <- y - cbind(1, dummies) %*% bdummies0
	  fit1 <- .robreg3S_noDummies(yt, x, filter=filter, alpha=alpha, ...)
	  b1 <- fit1$coef
	  fitdummies1 <- rlm( y - cbind(1, fit1$ximp) %*% b1 ~ dummies, method="M")
	  bdummies1 <- suppressWarnings(as.numeric(coef(fitdummies1)))
	  fit0 <- fit1
	  bdummies0 <- bdummies1
	}
	coeff.3S <- c(b1, bdummies1[-1])
	names(coeff.3S) <- c("(Intercept)", colnames(x), colnames(dummies))
	# yt <- y - dummies %*% bdummies1
	# acov.3S.nodummies <- .robreg3S.noDummies.acov(yt, fit0$ximp, coeff.3S[1:(1+p)], fit0$Syx, fit0$myx)
	# acov.3S.nodummies.se <- acov.3S.nodummies$se
	# acov.3S.nodummies <- acov.3S.nodummies$asv
	# acov.3S.dummies <- vcov(fitdummies1)*n
	# acov.3S <- matrix( NA, 1+p+pd, 1+p+pd)
	# acov.3S[1:(1+p), 1:(1+p)] <- acov.3S.nodummies
	# acov.3S[(2+p):(1+p+pd), (2+p):(1+p+pd)] <- acov.3S.dummies
	# se.3S <- sqrt(diag(acov.3S$asv)/n)
	# tab.3S.x <- cbind( coeff.3S[1:(1+p)], se.3S, coeff.3S[1:(1+p)]/se.3S, pnorm(abs(coeff.3S[1:(1+p)]/se.3S), lower.tail=FALSE)*2)
	acov.3S <- matrix( NA, 1+p+pd, 1+p+pd)
	acov.3S.nodummies.se <- NA
	tab.3S.x <- cbind(  coeff.3S[1:(1+p)], NA, NA, NA)
	tab.3S.dummies <- cbind( coeff.3S[colnames(dummies)], NA, NA, NA)
	tab.3S <- rbind(tab.3S.x, tab.3S.dummies) 
	colnames(tab.3S) <- c("Coef", "Asym.Std.Err.", "Z", "Pr(>|Z|)")
	ximp <- fit0$ximp
	colnames(ximp) <- colnames(x)
	result <- list(Summary.Table=tab.3S,
	               coef=coeff.3S,
	               acov=acov.3S,
	               resid=c(y - cbind(1, x) %*%b1 - dummies %*% bdummies1[-1]),
	               sigma.hat=acov.3S.nodummies.se,
	               MD=mahalanobis( cbind(y, x), fit0$myx, fit0$Syx),
	               weight=fit0$weight,
	               Syx=fit0$Syx,
	               myx=fit0$myx,
	               xfilter=xf,
	               ximpute=ximp
	)
	result
}




