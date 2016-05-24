factorScorePfa =
function(x, factors = 2, covmat = NULL, cor = FALSE, rotation = c("varimax", "none"), scoresMethod = c("none", "regression", "Bartlett")){
   cl = match.call()

   ## the default is computing things (factor scores etc) using the covariance matrix S!
   ## to use the correlation matrix, use cor = TRUE
   if (!is.null(covmat)) {
		if (is.list(covmat)) {
			S = covmat$cov # robust/classical covariance matrix
			center = covmat$center # robust/classical center
		}
		else { # covmat is not a list, e.g., covmat = S
			S = covmat # robust covariance matrix
			center = NULL # robust center
		}
	}
	else if (!is.null(x)) {
		S = cov(x) # classical covariance matrix
		center = colMeans(x) # classical center
	}
    else # covmat == NULL and x == NULL
		stop("no covmat or x provided")
	
	covariance = S
	correlation = cov2cor(S)
	
	if (cor == TRUE) {
		S = correlation # now S is the correlation matrix
	}

   d = 1/diag(solve(S))
   p = nrow(S); diag_S = diag(S); sum_rank = sum(diag_S)
   rowname = paste("X", 1:p, sep = "")
   colname = paste("Factor", 1:factors, sep = "")
   A0 = matrix(0, nrow = p, ncol = factors, 
             dimnames = list(rowname, colname))

   kmax = 20; k = 1; h = diag_S-d
   reducedCorrelation = S
   repeat{
      diag(reducedCorrelation) =  h; h1 = h; eig = eigen(reducedCorrelation)
      for (i in 1:factors)
         A0[,i] = sqrt(eig$values[i])*eig$vectors[,i]

      h = diag(A0 %*% t(A0))
      if ((sqrt(sum((h-h1)^2))<1e-4)|k == kmax) break
      k = k+1
   }

   if (missing(rotation) || rotation == "varimax")
        A = varimax(A0, normalize = TRUE)$loadings
   else if (rotation == "none")
        A = A0
   else cat("undefined rotation method, try rotation = 'varimax' or 'none' \n")
   
   # A is the factor loadings after rotation. The following for loop makes sure that 
   # in each column of A, the entry with the largest absolute value is always positive!
   for (i in 1:factors){
      if (A[,i][which.max(abs(A[,i]))]<0){
         A[,i] = -A[,i]
      }
   }

   h = diag(A%*%t(A))
   specific = diag_S-h

   if (missing(scoresMethod)) scoresMethod = "none"

   scoringCoef = F = meanF = corF = n.obs = NULL 
   if (!missing(x)) {
      n.obs = nrow(x)
	  scaledX = {if (cor == TRUE)
					scale(x, center = center, scale = sqrt(diag(covariance))) # standardized transformation, center and covariance maybe classical or robust
				 else # cor == FALSE
					scale(x, center = center, scale = FALSE) # centralized transformation, center maybe classical or robust
				}

      if (scoresMethod == "regression"){
		# compute the scoring coefficient
		scoringCoef = t(A) %*% solve(S)
		# compute scores
		F = scaledX %*% t(scoringCoef)
		# F = scale(x, center = center, scale = FALSE) %*% t(scoringCoef) # only do the centered transformation, to be compatible with the covariance matrix S.
		# F = scale(x, center = center, scale = sqrt(diag(S))) %*% t(scoringCoef) # standardized transformation using the robust center and scale, maybe incompatible with rrcov:::.myellipse()
		# F = scale(x) %*% t(scoringCoef) # standardized transformation
		# the sample mean of the scores F
		meanF = apply(F,2,mean)
		# the sample correlation matrix of the scores F
		corF = cor(F)
      }
   
      if (scoresMethod ==  "Bartlett"){
		# compute the scoring coefficient
		ADA.inv = solve(t(A) %*% diag(1/specific) %*% A)
		scoringCoef = ADA.inv %*% t(A) %*% diag(1/specific)
		# compute scores
		F = scaledX %*% t(scoringCoef)
		# F = scale(x, center = center, scale = FALSE) %*% t(scoringCoef) # only do the centered transformation, to be compatible with the covariance matrix S.
		# F = scale(x, center = center, scale = sqrt(diag(S))) %*% t(scoringCoef) # standardized transformation using the robust center and scale, maybe incompatible with rrcov:::.myellipse()
		# F = scale(x) %*% t(scoringCoef) # standardized transformation
		# the sample mean of the scores F
		meanF = apply(F,2,mean)
		# the sample correlation matrix of the scores F
		corF = cor(F)
      }
   }

   method = c("pfa") # Principal Factor Method
   res = list( call = cl,
		loadings = A, 
		communality = h,
		uniquenesses = specific,
		covariance = covariance, # robust/classical covariance matrix
		correlation = correlation, # robust/classical correlation matrix
		usedMatrix = S, # covariance or correlation matrix according to the value of cor
		reducedCorrelation = reducedCorrelation, # the last reduced correlation matrix
		factors = factors,
		method = method,
		scores = F,
		scoringCoef = scoringCoef,
		meanF = meanF, 
		corF = corF,
		scoresMethod = scoresMethod,
		n.obs = n.obs,
		center = center,
		eigenvalues = eigen(S)$values) # the eigenvalues of the usedMatrix S
   class(res) = "factorScorePfa"
   res
}
