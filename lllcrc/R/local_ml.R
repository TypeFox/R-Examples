#' Use odd-even formula to fit saturated LLM
#' 
#' See the odd-even formula in Fienberg 1972, or a localized version of it in
#' Zwane 2004.
#' 
#' The function is vectorized -- it gets applied to each row of the input
#' matrix, resulting in a vector of estimated rates of missingness or the total
#' number missing
#' 
#' @param dens Matrix with column names corresponding to the possible capture
#' patterns
#' @author Zach Kurtz
#' @references 
#' Fienberg SE (1972). "The Multiple Recapture Census for Closed
#' Populations and Incomplete $2^k$ Contingency Tables." \emph{Biometrika},
#' \bold{59}(3), pp. 591.
#' @export saturated.local
saturated.local = function(dens){ #dens = morph.ydens[i,]
 ## Given a data.frame of conditional densities, impute the missing cell
  Y = colnames(dens)
  eo = odd.evens(Y)
  evens = eo$even;  odds = eo$odd
  k = nchar(Y[1])
 # include some checks -- weird stuff should happen
 #   if any of the 2^k-1 expected terms are missing
  if(length(evens) < 2^(k-1)-1){out = 0/0	
 }else if(length(odds) < 2^(k-1)){out = 0
 }else{
  l.num = rowSums(log(dens[,odds, drop = FALSE]))
  l.den = rowSums(log(dens[,evens, drop = FALSE]))
  out = exp(l.num - l.den)}
  return(out)}



#' Determine the even-ness of capture patterns
#' 
#' Return a list containing the indices of capture patterns with an even number
#' of captures and the indices that of capture patterns with an odd number of
#' captures
#' 
#' 
#' @param Y A vector of capture patterns, such as the output of the function
#' \code{patterns.possible}
#' @author Zach Kurtz
#' @export odd.evens
odd.evens = function(Y){
 ## Given a vector of capture patterns strings, return the indices of
 #   the strings with an odd and even number of 1's, respectively
  sum.ci = function(k) sum(strsplit(Y[k], "")[[1]] == "1")
  odds.evens = sapply(1:length(Y), sum.ci)
  evens = which(odds.evens%%2 == 0)
  odds  = which(odds.evens%%2 != 0)
  eo = list(even = evens, odd = odds)
  return(eo)}


#' Maximum likelihood for log-linear coefficients
#' 
#' A simplified version of \code{glm} that does only parameter estimation
#' 
#' 
#' @param predictors The columns of the standard design matrix to include in
#' the model.  For example, "c1", "c2" for main effects, and "c12" for
#' interactions.
#' @param data A design matrix with cell counts included
#' @param normalized Logical: If TRUE, include a normalization step after
#' coefficient estimation, which resets the value of the intercept so that the
#' sum of predicted values is exactly 1
#' @param precision Controls the precision of the coefficient estimates.  A
#' higher number is less precise.  1 corresponds to machine epsilon.
#' @details Maximize the Poisson likelihood using BFGS in \code{optim()}.
#' @return The vector of estimated log-linear coefficients.  The first
#' coefficient is the intercept, and the remaining ones correspond to the
#' \code{predictors} argument, in that order
#' @author Zach Kurtz
#' @export zglm
zglm = function(predictors, data, normalized = TRUE, precision = 1000)
{  #predictors = best.terms, data = ddat, normalized = normalized
	# identify the relevant design matrix
	M = as.matrix(cbind(rep(1, nrow(data)), data[, predictors]))
	#initialize beta
	beta = rep(0, length(predictors)+1)
	#pick out the counts
	c = data[,"c"]
	#state the log-likelihood (up to a constant) for the data c and parameter beta:
	neg.pois.log.like.prop = function(beta){
		log.lambda = M%*%beta # log-expected cell counts under poisson model
		return(-sum(-exp(log.lambda) + c*log.lambda))
	}
	#state the gradient of the log-likelihood:
	grad.fun = function(beta){a = exp(M%*%beta)-c; return(t(a)%*%M)}
	#compute the MLE parameters
	beta = optim(beta, neg.pois.log.like.prop, method = "BFGS", gr = grad.fun, 
			control = list(reltol = precision*.Machine$double.eps))$par 
	if(normalized == TRUE)	beta[1] = beta[1] - log(sum(exp(M%*%beta)))
	out = beta # pred = exp(beta[1])
	return(out)
}

## A stripped-down high-speed version of Poisson regression with likelihood maximization by IRLS


#' Maximum likelihood for log-linear coefficients
#' 
#' A simplified version of \code{glm} that does only parameter estimation.
#' This attempts to mimic the IRLS routine invoked by \code{glm}, without
#' returning extra ``baggage" such as standard errors.
#' 
#' The main purpose of \code{pirls} is to obtain speed, for the special
#' circumstance in which one must fit gazillions of Poisson regression models,
#' where the only quantity of interested is the point estimates of the
#' regression coefficients.  Matrix inversion is one of the most time consuming
#' steps in the function, and the overall speed can be improved by about 20
#' percent by modifying the source code to replace the \code{solve()} command
#' with \code{.Internal(La_solve())}.
#' 
#' @param predictors The columns of the standard design matrix to include in
#' the model.  For example, "c1", "c2" for main effects, and "c12" for
#' interactions.
#' @param data A design matrix with cell counts included in the column named
#' "c".  Must be a matrix (not a data frame)!
#' @param epsilon Convergence tolerance, intended to play the same role as
#' \code{epsilon} the control parameters for \code{glm.fit}.
#' @param iter.max The maximum number of IRLS iterations.  A warning appears if
#' this maximum is ever reached.
#' @param normalized Logical: If TRUE, include a normalization step after
#' coefficient estimation, which resets the value of the intercept so that the
#' sum of predicted values is exactly 1
#' @return The vector of estimated log-linear coefficients.  The first
#' coefficient is the intercept, and the remaining ones correspond to the
#' \code{predictors} argument, in that order
#' @author Zach Kurtz
#' @export pirls
pirls = function(predictors, data, epsilon = 1e-8, iter.max = 25, normalized = TRUE) 
{
	p = length(predictors) + 1
	y = data[,"c"]
	M = cbind(rep(1, nrow(data)), data[, predictors])
	beta = c(log(mean(y)), rep(0, p-1))
	Mbeta = M%*%beta
	mu = exp(Mbeta)
	dev = -2*(y%*%Mbeta  - sum(mu))
	iter = 0
	crit = 1
	while((iter < iter.max) & (epsilon < crit)){
		iter = iter+1
		W = diag(as.numeric(mu))
		tM = t(M)
		a = (tM%*%W) %*% M
		b = diag(1, nrow(a)) #only a device for setting up La_solve to get the inverse
		#Slightly faster:  inv = .Internal(La_solve(a, b, tol = .Machine$double.eps))
		inv = solve(a, b, tol = .Machine$double.eps)
		beta = beta + inv %*% (tM%*%(y - mu))  #irls iterative update
		Mbeta = M%*%beta
		mu = exp(Mbeta)
		newdev = -2*(y%*%Mbeta  - sum(mu))
		crit = abs(newdev-dev)/(abs(dev) + 0.1)
		dev = newdev
	}
	if(iter == iter.max) warning("pirls failed to converge within 25 iterations")
	if(normalized == TRUE) beta[1] = beta[1] - log(sum(exp(Mbeta)))
	return(beta)
}


#' Maximum likelihood estimation for fixed LLLMs
#' 
#' A workhorse function for \code{apply.local.ml}.  For a fixed log-linear model,
#' get best fit using maximum local likelihood estimation.
#' 
#' Specify exactly one of the two arguments \code{predictors}, \code{ml.meth}
#' 
#' @param pdat A design matrix with cell counts (possibly fractional) 
#' included in the column named "c".
#' 
#' @param ml.meth A model specification such as "indep" or "indep.mono".
#' @param predictors A character vector of predictors of the form "c1", "c2"
#' for main effects, or "c12" for an interaction.
#' @param k Number of lists. 
#' @author Zach Kurtz
#' @references 
#' Fienberg SE (1972). "The Multiple Recapture Census for Closed
#' Populations and Incomplete $2^k$ Contingency Tables." \emph{Biometrika},
#' \bold{59}(3), pp. 591.
#' @references 
#' Cormack RM (1989). "Log-linear models for capture-recapture."
#' \emph{Biometrics}, pp. 395-413.
#' @export local.ml
local.ml = function(pdat, ml.meth = NULL, predictors = NULL, k){ 
	if(is.null(ml.meth)&is.null(predictors)){
		stop("Missing both ml.meth and predictors - specify exactly one of these")
	}else
	if(!is.null(ml.meth)&!is.null(predictors)){
		stop("There is both ml.meth and predictors - specify only one of these")
	}else
	# specify any formula with fields in pdat
	if(!is.null(predictors)){
		beta = pirls(predictors = predictors, data = as.matrix(pdat), normalized = FALSE)
	}else 
	if(ml.meth == "indep"){# general independence model
		terms = paste0(rep("c",k), 1:k) #paste0(names(pdat)[1:k], collapse = "+")
		beta = pirls(predictors = terms, data = as.matrix(pdat), normalized = FALSE)
	}else  # independence special case: no list effects
	if(ml.meth == "indep.mono"){
	pdat$list.count = rowSums(pdat[,1:k])
	beta = pirls(predictors = "list.count", data = as.matrix(pdat), normalized = FALSE) 
	}else{stop("ml.meth given is not yet implemented in ml.local")
	}
	pred = exp(beta[1])
	return(pred)
}



#' Fit LLLMs
#' 
#' Fit a previously selected log-linear model at each covariate vector.
#' 
#' This function is closely related to \code{apply.ic.fit}, with the difference
#' being that the user supplies the log-linear model to be applied locally
#' (instead of using an information criterion to select a potentially different
#' model at each point).
#' 
#' @param dens Same as \code{ydens} in apply.ic.fit.
#' @param ml.meth A model specification, in case \code{predictors} is NOT
#' specified.  Options include "indep" and "indep.mono".
#' @param predictors A character vector of predictors of the form "c1", "c2"
#' for main effects, or "c12" for an interaction.
#' @param design.style The kind of model design; Rasch or not Rasch
#' @return A vector of estimated rates of missingness corresponding to each
#' point.
#' @author Zach Kurtz
#' @export apply.local.ml
apply.local.ml = function(dens, ml.meth = NULL, predictors = NULL, design.style = "standard"){  #dens = tr$ydens
	# dens = dens; ml.meth = "formulas"; formula = "c~c12+c13+c23"
	dens = data.matrix(dens) #if dens is a data.frame, converts to matrix
	k = nchar(colnames(dens)[1])
	if(design.style == "standard") design.array = make.design.matrix(k)
	if(design.style == "rasch") design.array = make.design.matrix(k, rasch = TRUE)
	if(!is.element(design.style, c("standard", "rasch"))) stop("invalid design.style in apply.local.ml")
	get.row.i = function(i) { #i = 10
		design.array$c = as.numeric(dens[i,])
		out = local.ml(design.array, ml.meth = ml.meth, predictors = predictors, k)
		return(out)
	}
	res = sapply(1:nrow(dens), get.row.i)
	return(res)
}
