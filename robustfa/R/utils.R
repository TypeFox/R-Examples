detail = function(x){
res  =  list(x = x,
	    isS4 = isS4(x),
	    isObject = is.object(x),
	    class = class(x),
	    attributes = attributes(x))
res
}

computeScores = function(out, x = data, covmat = covmat, cor = cor, scoresMethod = scoresMethod) {
    if (!is.list(out))
	stop("out is not a list")
	
	## the default is computing things (factor scores etc) using the covariance matrix S!
	## to use the correlation matrix, use cor = TRUE
	if (!is.null(covmat)) {
		if (is.list(covmat)) {
			S = covmat$cov # robust/classical covariance matrix
			center = covmat$center # robust/classical center
		}
		else { # covmat is not a list, e.g., covmat = S
			S = covmat # robust/classical covariance matrix
			center = NULL # robust/classical center
		}
	}
	
	covariance = S
	correlation = cov2cor(S)
	
	if (cor == TRUE) {
		S = correlation # now S is the correlation matrix
	}
	
	scaledX = {if (cor == TRUE)
					scale(x, center = center, scale = sqrt(diag(covariance))) # standardized transformation, center and covariance maybe classical or robust
				 else # cor == FALSE
					scale(x, center = center, scale = FALSE) # centralized transformation, center maybe classical or robust
				}

    if (scoresMethod == "none"){
          scoringCoef = F = meanF = corF = NULL
    }
    else if (scoresMethod == "regression"){
		# compute the scoring coefficient
		scoringCoef  = t(out$loadings[]) %*% solve(S)
		# compute scores
		F = scaledX %*% t(scoringCoef)
		# F = scale(x, center = center, scale = FALSE) %*% t(scoringCoef) # only do the centered transformation, to be compatible with the covariance matrix S.
		# F = scale(x, center = center, scale = sqrt(diag(S))) %*% t(scoringCoef) # standardized transformation using the robust center and scale, maybe incompatible with rrcov:::.myellipse()
		# F = scale(x) %*% t(scoringCoef) # standardized transformation
		# the sample mean of the scores F
		meanF  = apply(F, 2, mean)
		# the sample correlation matrix of the scores F
		corF  = cor(F)
	}
	else{ ## (scoresMethod == "Bartlett")
		# compute the scoring coefficient
		ADA.inv  = solve(t(out$loadings[]) %*% diag(1/out$uniquenesses) %*% out$loadings[])
		scoringCoef  = ADA.inv %*% t(out$loadings[]) %*% diag(1/out$uniquenesses)
		# compute scores
		F = scaledX %*% t(scoringCoef)
		# F = scale(x, center = center, scale = FALSE) %*% t(scoringCoef) # only do the centered transformation, to be compatible with the covariance matrix S.
		# F = scale(x, center = center, scale = sqrt(diag(S))) %*% t(scoringCoef) # standardized transformation using the robust center and scale, maybe incompatible with rrcov:::.myellipse()
		# F = scale(x) %*% t(scoringCoef) # standardized transformation
		# the sample mean of the scores F
		meanF  = apply(F, 2, mean)
		# the sample correlation matrix of the scores F
		corF  = cor(F)
	}
	 
	res = out
	res$scoringCoef = scoringCoef
	res$scores = F
	res$meanF = meanF
	res$corF = corF
	res$eigenvalues = eigen(S)$values # the eigenvalues of the usedMatrix S
	res$covariance = covariance
	res$correlation = correlation
	res$usedMatrix = S
	res$reducedCorrelation = NULL
	res
}

compute_cov_cor = function(x, control){

cov_x = CovRobust(x = x, control = control)
cov_scale_x = CovRobust(x = scale(x), control = control)

## S_r != S_r_tilda? Yes!
## R_r == R_r_tilda? Yes!
res = list(
		S_r = getCov(cov_x),
		S_r_tilda = getCov(cov_scale_x),
		R_r = getCorr(cov_x),
		R_r_tilda = getCorr(cov_scale_x))
res
}

