print.summary.aodml <- function(x, ..., digits = max(3, getOption("digits") - 3)) {
	
	Object <- x$object
	cat("\nCall:\n")
	print(Object$call)

	# checks whether convergence problems occurred
  n <- Object$code
  iter <- Object$iterations
  msg <- Object$msg
  if(n == 0)
    cat("\nConvergence was obtained after " , iter, " iterations.\n", sep = "")

	## Vector b
	
	# print estimated coef, if any
  BCoef <- x$BCoef
  if(nrow(BCoef) > 0){
    nam <- rownames(BCoef)
  	List <- vector(mode = "list", length = 4)
  	for(i in 1:4){
      X <- BCoef[, i]
## Modif RL 16/06/2013
##      List[[i]] <- format(X, digits = digits, scientific = FALSE)
      List[[i]] <- sapply(X, function(x){
        sci <- if(abs(x) < 10^(-4)) TRUE else FALSE
        format(x, nsmall = 4, digits = digits, scientific = sci)
        })
  	  }
    BCoeftext <- as.data.frame(t(do.call("rbind", List)))
    rownames(BCoeftext) <- nam
    colnames(BCoeftext) <- c("Estimate", "Std. Error", "z value", "Pr(> |z|)")
    cat("\nMu coefficients:\n")
    print(BCoeftext)
  }

	# print coef set to a fixed value, if any
  FixedBCoef <- x$FixedBCoef
  if(nrow(FixedBCoef) > 0) {
    cat("\nMu coefficients (b) set to fixed values:\n")
    print(FixedBCoef)
  }

	## Vector phi
	# compute new var-cov mat, coef vector and position of term(s) to be tested
  Phi <- x$Phi
  if(nrow(Phi) > 0) {
    nam <- rownames(Phi)
    List <- vector(mode = "list", length = 2)
  	for(i in 1:2) {
  		X <- Phi[,i]
  		## Modif RL 16/06/2013
  		##      List[[i]] <- format(X, digits = digits, scientific = FALSE)
  		List[[i]] <- sapply(X, function(x){
  		  sci <- if(abs(x) < 10^(-4)) TRUE else FALSE
  		  format(x, nsmall = 4, digits = digits, scientific = sci)
  		})
  	}
  	Phitext <- as.data.frame(t(do.call("rbind", List)))
  	rownames(Phitext) <- nam
  	colnames(Phitext) <- c("Estimate", "Std. Error") #, "z value", "Pr(> z)"
  	cat("\nPhi coefficients: (scale = ",  Object$phi.scale, ")\n", sep = "")
  	print(Phitext)
  }

	# print coef set to a fixed value, if any
  FixedPhi <- x$FixedPhi
  if(nrow(FixedPhi) > 0){
  	cat("\nPhi coefficients set to fixed values:\n")
  	print(FixedPhi)
  }
  akic <- AIC(Object)
  ll    <- format(Object$logL, digits = digits, scientific = FALSE)
  dfm <- format(Object$df.model)
  dfres <- format(df.residual(Object))
  dev   <- format(deviance(Object), digits = digits, scientific = FALSE)
  aic   <- format(akic[, 3], digits = digits, scientific = FALSE)
  aicc  <- format(akic[, 4], digits = digits, scientific = FALSE)
  res   <- c(ll, dfm, dfres, dev, aic, aicc)
  names(res) <- c("Log-lik", "df.model", "df.resid", "Deviance", "AIC", "AICc")
  cat("\nLog-likelihood statistics\n")
  print(res, quote = FALSE)
  invisible(x)
  }
