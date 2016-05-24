#' Fit distributions to elicited probabilities 
#' 
#' Takes elicited probabilities as inputs, and fits parametric distributions
#' using least squares on the cumulative distribution function. If separate
#' judgements from multiple experts are specified, the function will fit one
#' set of distributions per expert.
#' 
#' 
#' @param vals A vector of elicited values for one expert, or a matrix of
#' elicited values for multiple experts (one column per expert). Note that the
#' an elicited judgement about X should be of the form P(X<= vals[i,j]) =
#' probs[i,j]
#' @param probs A vector of elicited probabilies for one expert, or a matrix of
#' elicited values for multiple experts (one column per expert). A single
#' vector can be used if the probabilities are the same for each expert. For
#' each expert, the smallest elicited probability must be less than 0.4, and
#' the largest elicited probability must be greater than 0.6.
#' @param lower A single lower limit for the uncertain quantity X, or a vector
#' of different lower limits for each expert. Specifying a lower limit will
#' allow the fitting of distributions bounded below.
#' @param upper A single upper limit for the uncertain quantity X, or a vector
#' of different lower limits for each expert. Specifying both a lower limit and
#' an upper limit will allow the fitting of a Beta distribution.
#' @param weights A vector or matrix of weights corresponding to vals if
#' weighted least squares is to be used in the parameter fitting.
#' @param tdf The number of degrees of freedom to be used when fitting a
#' t-distribution.
#' @return 
#' \item{Normal}{Parameters of the fitted normal distributions.}
#' \item{Student.t}{Parameters of the fitted t distributions. Note that (X -
#' location) / scale has a standard t distribution. The degrees of freedom is
#' not fitted; it is specified as an argument to \code{fitdist}.}
#' \item{Gamma}{Parameters of the fitted gamma distributions. Note that E(X) =
#' shape / rate.} 
#' \item{Log.normal}{Parameters of the fitted log normal
#' distributions: the mean and standard deviation of log X.}
#' \item{Log.Student.t}{Parameters of the fitted log student t distributions.
#' Note that (log(X) - location) / scale has a standard t distribution. The
#' degrees of freedom is not fitted; it is specified as an argument to
#' \code{fitdist}.} 
#' \item{Beta}{Parameters of the fitted beta distributions. X
#' is scaled to the interval [0,1] via Y = (X - \code{lower})/(\code{upper} -
#' \code{lower}), and E(Y) = shape1 / (shape1 + shape2).} 
#' \item{ssq}{Sum of
#' squared errors for each fitted distribution and expert. Each error is the
#' different between an elicited cumulative probability and the corresponding
#' fitted cumulative probability.} 
#' \item{best.fitting}{The best fitting
#' distribution for each expert, determined by the smallest sum of squared
#' errors.} 
#' \item{vals}{The elicited values used to fit the distributions.}
#' \item{probs}{The elicited probabilities used to fit the distributions.}
#' \item{limits}{The lower and upper limits specified by each expert (+/- Inf
#' if not specified).}
#' @note The least squares parameter values are found numerically using the
#' \code{optim} command. Starting values for the distribution parameters are
#' chosen based on a simple normal approximation: linear interpolation is used
#' to estimate the 0.4, 0.5 and 0.6 quantiles, and starting parameter values
#' are chosen by setting E(X) equal to the 0.5th quantile, and Var(X) = (0.6
#' quantile - 0.4 quantile)^2 / 0.25. Note that the arguments \code{lower} and
#' \code{upper} are not included as elicited values on the cumulative
#' distribution function. To include a judgement such as P(X<=a)=0, the values
#' a and 0 must be included in \code{vals} and \code{probs} respectively.
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @examples
#' \dontrun{
#' # One expert, with elicited probabilities
#' # P(X<20)=0.25, P(X<30)=0.5, P(X<50)=0.75
#' # and X>0.
#' v <- c(20,30,50)
#' p <- c(0.25,0.5,0.75)
#' fitdist(vals=v, probs=p, lower=0)
#' 
#' # Now add a second expert, with elicited probabilities
#' # P(X<55)=0.25, P(X<60=0.5), P(X<70)=0.75
#' v <- matrix(c(20,30,50,55,60,70),3,2)
#' p <- c(0.25,0.5,0.75)
#' fitdist(vals=v, probs=p, lower=0)
#' 
#' # Two experts, different elicited quantiles and limits.
#' # Expert 1: P(X<50)=0.25, P(X<60=0.5), P(X<65)=0.75, and provides bounds 10<X<100
#' # Expert 2: P(X<40)=0.33, P(X<50=0.5), P(X<60)=0.66, and provides bounds 0<X
#' v <- matrix(c(50,60,65,40,50,60),3,2)
#' p <- matrix(c(.25,.5,.75,.33,.5,.66),3,2)
#' l <- c(10,0)
#' u <- c(100, Inf)
#' fitdist(vals=v, probs=p, lower=l, upper=u)
#' }
#' @import stats
#' @export
#' 
#' 
fitdist <-
function(vals, probs, lower = -Inf, upper = Inf, weights = 1, tdf = 3){
	
	if(is.matrix(vals)==F){vals<-matrix(vals, nrow = length(vals), ncol = 1)}
	if(is.matrix(probs)==F){probs <- matrix(probs, nrow = nrow(vals), ncol = ncol(vals))}
	if(is.matrix(weights)==F){weights <- matrix(weights, nrow = nrow(vals), ncol = ncol(vals))}
	if(length(lower)==1){lower <- rep(lower, ncol(vals))}
	if(length(upper)==1){upper <- rep(upper, ncol(vals))}
	if(length(tdf)==1){tdf <- rep(tdf, ncol(vals))}
	
  
	
	n.experts <- ncol(vals)
	normal.parameters <- matrix(NA, n.experts, 2)
	t.parameters <- matrix(NA, n.experts, 3)
	gamma.parameters <- matrix(NA, n.experts, 2)
	lognormal.parameters <- matrix(NA, n.experts, 2)
	logt.parameters <- matrix(NA, n.experts, 3)
	beta.parameters <- matrix(NA, n.experts, 2)
	ssq<-matrix(NA, n.experts, 6)
	expertnames <- paste("expert.", 1:n.experts, sep="")
	
	limits <- data.frame(lower = lower, upper = upper)
	row.names(limits) <- expertnames
	
	for(i in 1:n.experts){
		if (min(probs[,i]) > 0.4 ){stop("smallest elicited probability must be less than 0.4")}
		if (min(probs[,i]) < 0 | max(probs[,i]) > 1 ){stop("probabilities must be between 0 and 1")}
		if (max(probs[,i]) < 0.6 ){stop("largest elicited probability must be greater than 0.6")}
    if (min(vals[,i]) < lower[i]){stop("elicited parameter values cannot be smaller than lower parameter limit")}
		if (max(vals[,i]) > upper[i]){stop("elicited parameter values cannot be greater than upper parameter limit")}
		if (tdf[i] <= 0 ){stop("Student-t degrees of freedom must be greater than 0")}
		if (min(probs[-1,i] - probs[-nrow(probs),i]) <= 0 ){stop("probabilities must be specified in ascending order")}
		if (min(vals[-1,i] - vals[-nrow(vals),i]) <= 0 ){stop("parameter values must be specified in ascending order")}
    
	  minprob <- min(probs[, i])
	  maxprob <- max(probs[, i])
		
	  q.fit <- approx(x = probs[,i], y = vals[,i], xout = c(0.4, 0.5, 0.6))$y
	  l <- q.fit[1]
	  u <- q.fit[3]
	  
	  if(minprob > 0 & maxprob < 1){
		  minvals <- min(vals[, i])
		  maxvals <- max(vals[, i])
		  minq <- qnorm(minprob)
		  maxq <- qnorm(maxprob)
		  m <- (minvals * maxq - maxvals * minq) / (maxq - minq)
		  v <- ((maxvals - minvals) / (maxq - minq))^2
		}else{
		  m <- q.fit[2]
		  v<- (u - l)^2 / 0.25
		} 
	
		normal.fit <- optim(c(m, 0.5*log(v)), normal.error, values = vals[,i], probabilities = probs[,i], weights = weights[,i])   
    normal.parameters[i,] <- c(normal.fit$par[1],exp(normal.fit$par[2]))
    ssq[i,1] <- normal.fit$value
	
		t.fit <- optim(c(m, 0.5*log(v)), t.error, values = vals[,i], probabilities = probs[,i], weights = weights[,i], degreesfreedom = tdf[i])
    	t.parameters[i, 1:2] <- c(t.fit$par[1], exp(t.fit$par[2]))
    	t.parameters[i, 3] <- tdf[i]
    	ssq[i,2] <- t.fit$value
	
	
		if(lower[i] > -Inf){
			vals.scaled1 <- vals[,i] - lower[i]
			m.scaled1 <- m - lower[i]
		
			gamma.fit<-optim(c(log(m.scaled1^2/v), log(m.scaled1/v)), gamma.error, values = vals.scaled1, probabilities = probs[,i], weights = weights[,i])
    		gamma.parameters[i,] <- exp(gamma.fit$par)
    		ssq[i,3] <- gamma.fit$value
    		
    		std<-((log(u)-log(l))/1.35)
    	
    		lognormal.fit <- optim(c(log(m.scaled1), log(std)), lognormal.error, values = vals.scaled1, probabilities = probs[,i], weights = weights[,i])
    		lognormal.parameters[i, 1:2] <- c(lognormal.fit$par[1], exp(lognormal.fit$par[2]))
    		ssq[i,4] <- lognormal.fit$value
    	
    		logt.fit <- optim(c(log(m.scaled1), log(std)), logt.error, values = vals.scaled1, probabilities = probs[,i], weights = weights[,i], degreesfreedom = tdf[i])
    		logt.parameters[i,1:2] <- c(logt.fit$par[1], exp(logt.fit$par[2]))
    		logt.parameters[i,3] <- tdf[i]
    		ssq[i,5] <- logt.fit$value
		}
	
		if((lower[i] > -Inf) & (upper[i] < Inf)){
			vals.scaled2 <- (vals[,i] - lower[i]) / (upper[i] - lower[i])
			m.scaled2 <- (m - lower[i]) / (upper[i] - lower[i])
			v.scaled2 <- v / (upper[i] - lower[i])^2
		
			alp <- abs(m.scaled2 ^3 / v.scaled2 * (1/m.scaled2-1) - m.scaled2)
    		bet <- abs(alp/m.scaled2 - alp)
    		beta.fit <- optim(c(log(alp), log(bet)), beta.error, values = vals.scaled2, probabilities = probs[,i], weights = weights[,i])
    		beta.parameters[i,] <- exp(beta.fit$par)
    		ssq[i,6] <- beta.fit$value	
		
		}
	}
	dfn <- data.frame(normal.parameters)
	names(dfn) <-c ("mean", "sd")
	row.names(dfn) <- expertnames
		
	dft <- data.frame(t.parameters)
	names(dft) <-c ("location", "scale", "df")
	row.names(dft) <- expertnames
	
	dfg <- data.frame(gamma.parameters)
	names(dfg) <-c ("shape", "rate")
	row.names(dfg) <- expertnames
	
	dfln <- data.frame(lognormal.parameters)
	names(dfln) <-c ("mean.log.X", "sd.log.X")
	row.names(dfln) <- expertnames
	
	dflt <- data.frame(logt.parameters)
	names(dflt) <-c ("location.log.X", "scale.log.X", "df.log.X")
	row.names(dflt) <- expertnames
	
	dfb <- data.frame(beta.parameters)
	names(dfb) <-c ("shape1", "shape2")
	row.names(dfb) <- expertnames
	
	ssq <- data.frame(ssq)
	names(ssq) <- c("Normal", "Student-t", "Gamma", "Log normal", "Log Student-t", "Beta")
	row.names(ssq) <- expertnames
	
	index <- apply(ssq, 1, which.min)
	best.fitting <- data.frame(best.fit=names(ssq)[index])
	row.names(best.fitting) <- expertnames
  
	vals <- data.frame(vals)
	names(vals) <- expertnames
	
	probs <- data.frame(probs)
	names(probs) <- expertnames
		
  list(Normal = dfn, Student.t = dft, Gamma = dfg, Log.normal = dfln, Log.Student.t = dflt, Beta = dfb, ssq = ssq, best.fitting = best.fitting, vals = t(vals), probs = t(probs), limits = limits)
}
