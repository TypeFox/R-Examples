## Implement stepwise regression for local model selection by some Information Criterion (ic)


#' Compute an information criterion
#' 
#' Given the fitted parameter values for a log-linear model, compute an
#' information criterion.
#' 
#' Computes the conditional multinomial likelihood and uses it to compute the
#' specified information criterion
#' 
#' @param predictors A character vector of predictors of the form "c1", "c2"
#' for main effects, or "c12" for an interaction. The predictors to be used in
#' a log-linear model.  For example, "c1", "c2" for main effects, or "c12" for
#' an interaction.
#' @param ddat A data frame that is the design matrix for a log-linear model.
#' @param ic The information criterion to be computed.  Currently the AIC,
#' AICc, BIC, BICpi are implemented.
#' @param beta The vector of log-linear coefficients that were previously
#' estimated.
#' @return The value of the information criterion
#' @references Thesis of Zach Kurtz (2014), Carnegie Mellon University, Statistics
#' @references
#' Anderson DR and Burnham KP (1999). "Understanding information criteria
#' for selection among capture-recapture or ring recovery models." \emph{Bird
#' Study}, \bold{46}(S1), pp. S14-S21.
#' @author Zach Kurtz
#' @export get.IC
get.IC = function(predictors, ddat, ic, beta){ #mod = mod.pmml
	# Given data and a model, compute the IC.
	# Compute the probability array implied by the model:
	n = sum(ddat$c) 
	M = as.matrix(cbind(rep(1, nrow(ddat)), ddat[, predictors]))
	Ec = exp(M%*%beta)
	log.p = log(Ec/sum(Ec))
	q = length(predictors) # + 1 #DO NOT add in the intercept, as it is not a free parameter
	# This part makes sense if you carefully write out the likelihood function
	logL.choose.section = lfactorial(n)- sum(lfactorial(ddat$c))
	logL.kernel = sum(ddat$c*log.p)
	logL = logL.choose.section + logL.kernel

	if(ic == "AICc"){ out = -2*logL + q*2 + 2*q*(q+1)/(n - q - 1)
	}else if(ic == "BICpi"){ out = -2*logL + q*log(n/(2*pi))
	}else if(ic == "BIC"){ out = -2*logL + q*log(n)
	}else if(ic == "AIC"){ out = -2*logL + q*2
	}else{stop("Define your ic! -- says get.IC()")}
	return(out)
}



#' Compute an IC for several LLMs
#' 
#' Given several sets of log-linear terms, compute the IC corresponding to each
#' model.
#' 
#' 
#' @param models A list of character vectors, with each vector containing
#' column names from the associated log-linear design matrix.
#' For example, see the output of \code{\link{make.hierarchical.term.sets}()}.
#' @param ddat The log-linear design matrix.
#' @param ic The information criterion, such as AIC, AICc, BIC, or BICpi.
#' @param normalized Logical: TRUE means that beta0 will be adjusted so that
#' the log-linear model corresponds to cell probabilities instead of expected
#' cell counts.
#' @return A matrix with as many rows as there are entries in \code{models}.
#' The columns contain the point estimates of the population size, the
#' information criterion scores, and the information criterion weights for all
#' the models, which sum to one
#' @author Zach Kurtz
#' @export ic.all
ic.all = function(models, ddat, ic, normalized = normalized){ #mod = mod0; 
	# Use IC to select a local log-linear model from one of the entries in term.sets
	n.mod = length(models)
	results = matrix(NA, nrow = n.mod, ncol = 3)
	colnames(results) = c("est", "score", "wghts")
	mdat = as.matrix(ddat)
	for(i in 1:n.mod){ 
		predictors = models[[i]]
		beta = pirls(predictors, mdat, normalized = normalized)
		results[i,"score"] = get.IC(models[[i]], ddat, ic, beta)
		results[i,"est"] = exp(beta[1])
	}
	results[,"wghts"] = ic.wghts(results[,"score"])
	return(results)
} 



#' Select and fit an LLM
#' 
#' Use an information criterion to select a local log-linear model
#' 
#' Just like \code{flat.IC} except that it is designed to take in a local
#' average instead of a full capture-recapture dataset
#' 
#' @param densi A matrix with one row and 2^k-1 column containing cell counts
#' or empirical cell probabilities corresponding to all the possible capture
#' patterns.
#' @param models A list of character vectors, with each vector containing
#' column names from the associated log-linear design matrix.
#' For example, see the output of \code{\link{make.hierarchical.term.sets}()}.
#' @param N If you multiply \code{densi} by \code{N} and then sum over the
#' resulting vector, you should get the effective sample size.
#' @param ic The information criterion, such as AIC, AICc, BIC, or BICpi.
#' @param averaging Logical: TRUE means that we use information criterion
#' scores to do model averaging.
#' @param normalized Logical: TRUE means that beta0 will be adjusted so that
#' the log-linear model corresponds to cell probabilities instead of expected
#' cell counts.
#' @param rasch Logical: TRUE means that the Rasch model is a candidate.
#' @return \item{pred}{Estimated rate of missingness for the selected model}
#' \item{form}{Formula of the selected model}
#' @author Zach Kurtz
#' @export ic.fit
ic.fit = function(densi, models, N, ic, averaging = averaging, normalized = TRUE, rasch = FALSE){ 
	#i = 2; N = mdf[i]; densi = ydens[i,]; order.max = 2; stepwise = FALSE
	# Use an IC to select a local log-linear model. 
	ddat = string.to.array(densi, rasch = rasch) #ddat stands for "design data"
	ddat$c = ddat$c*N # If the data is 
	# Perform model selection using the ic, optionally stepwise
	icd = ic.all(models, ddat = ddat, ic, normalized = normalized)
	winner = which.min(icd[,"score"])
	best.terms = models[[winner]]
 	pred = ifelse(averaging, sum(icd[,"est"]*icd[,"wghts"]), icd[winner,"est"])
	# Return the winning formula and an estimate
	out = list(pred = pred, form = paste(best.terms, collapse = "+"))
	return(out)}



#' Select an LLLM at each point
#' 
#' Select LLLMs for each row of the input data.
#' 
#' See Kurtz (2013).  Each row of \code{ydens} corresponds to a covariate
#' vector, and contains a local average of multinomial capture pattern outcomes
#' across nearby points.  \code{apply.ic.fit} applies the function
#' \code{ic.fit} at each row.  The vector of local effective sample sizes is
#' crucial, and is specified in the \code{ess} argument.
#' 
#' @param ydens A matrix with 2^k-1 columns, one for each capture pattern.
#' Each row sums to 1; these are empirical capture pattern probabilities.
#' @param models A list of character vectors, with each vector containing
#' column names from the associated log-linear design matrix.
#' For example, see the output of \code{\link{make.hierarchical.term.sets}()}.
#' @param ess A vector of effective sample sizes, one for each row of ydens.
#' @param mct The number of population units that were observed for each row of
#' ydens.
#' @param ic The chosen information criterion.  Currently implemented: "AIC",
#' "AICc", "BIC", "BICpi".
#' @param cell.adj Logical: TRUE means that the cell adjustment of Evans and
#' Bonet (1995) is applied.
#' @param averaging Logical: TRUE means that the information criterion weights
#' are used to do model averaging, locally.
#' @param loud Logical: TRUE means that the progress is noted by printing the
#' number of the row of ydens currently being processed.
#' @return \item{lll}{An object of class "lllcrc"}
#' @author Zach Kurtz
#' @references Kurtz (2013)
#' @export apply.ic.fit
apply.ic.fit = function(ydens, models, ess, mct, ic, cell.adj, averaging, loud = TRUE){ 
 #ydens=sdat$hpi; ess=sdat$ess[,1]; mct = cdt[,"mct"]; ic = "AICc"; cell.adj = TRUE; averaging = FALSE; loud = TRUE
	if(cell.adj){
		k=nchar(colnames(ydens)[1])
		ydens = ydens + (1/2^(k-1))/ess
		ydens = ydens/rowSums(ydens)
	}
	ith.row = function(i) {
		if(loud) print(i)
		densi = ydens[i,,drop = FALSE] 
		out = ic.fit(densi, models, N = ess[i], ic, averaging = averaging)
		return(out)
	}
	out = t(sapply(1:nrow(ydens), ith.row))
	out = data.frame(out)
	out$form = unlist(out$form)
	out$pred = unlist(out$pred)
	loc.pred = out$pred*mct
	est = sum(loc.pred)
	loc.pred[mct == 0]=NA
	res = list(est = est, cpi0 = loc.pred, llform = out$form)
	return(res)
}

# Compute IC weights (i.e., AICc weights, or BIC, whatever)


#' Information criterion model weights
#' 
#' Compute model weights according to the information criterion scores of each
#' model.
#' 
#' The formula is quite simple: Identify the smallest (best) score among the
#' various models.  Subtract this minimum value from all of the scores, and
#' call the resulting set of scores $s$.  Compute exp(-0.5 s) for all the
#' scores, and normalize the resulting vector to obtain the vector of model
#' weights
#' 
#' @param scores The information criterion scores.
#' @return A vector of weights, which can be interpreted (loosely) as the
#' relative desireability of the models corresponding to the weights
#' @author Zach Kurtz
#' @references Burnham and Anderson (2002)
#' @export ic.wghts
ic.wghts = function(scores) {	
	scores = scores - min(scores)
	expsc = exp(-0.5*scores)
	denom = sum(expsc)
	wghts = expsc/denom
	return(wghts)}

