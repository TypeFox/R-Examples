

#' Compute estimates of the sampling weights of the respondent's observations
#' based on various estimators
#' 
#' 
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param weight.type A string giving the type of estimator to use. The options
#' are \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I/DS"},
#' and \code{"Arithemic Mean"}. It defaults to \code{"Gile's
#' SS"}.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{population.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing, the weights will sum to 1. Note that
#' this parameter is required for Gile's SS.
#' @param subset A logical expression subsetting rds.data.
#' @param control A list of control parameters for algorithm
#' tuning. Constructed using\cr
#' \code{\link{control.rds.estimates}}.
#' @param ... Additional parameters passed to the individual weighting algorithms.
#' @return A vector of weights for each of the respondents. It is of the same
#' size as the number of rows in \code{rds.data}.
#' @seealso \code{\link{rds.I.weights}}, \code{\link{gile.ss.weights}}, \code{\link{vh.weights}}
#' @export
compute.weights <- function(rds.data,
		weight.type = c("Gile's SS","RDS-I","RDS-I (DS)","RDS-II","Arithmetic Mean","Good-Fellows"),
		N = NULL, subset=NULL, control=control.rds.estimates(), ...){
	if(!is.null(N) && N < nrow(rds.data)){
		stop(sprintf("The population size, %d, is less than the sample
								size, %d. The population size must be at least as large as the
								sample size for the estimate to make sense.",N,nrow(rds.data)))
	}
	
	n <- nrow(rds.data)
	if(is.null(attr(rds.data,"network.size.variable")))
		stop("rds.data must have a network.size attribute.") #wait-is network.size a vector?
	init.deg <- as.numeric(rds.data[[attr(rds.data,"network.size.variable")]])
	deg = init.deg
	remvalues <- deg==0 | is.na(deg)
	if(any(remvalues)){
		warning(paste(sum(remvalues),"of",nrow(rds.data),
						"network sizes were missing or zero. The estimator will presume these are",max(deg,na.rm=TRUE)), call. = FALSE)
		deg[remvalues] <- max(deg,na.rm=TRUE)
	}
	if(missing(weight.type)){
		weight.type <- "Gile's SS"
	}
	weight.type <- match.arg(weight.type,
			c("Gile's SS","RDS-I", "RDS-II", "RDS-I (DS)","Arithmetic Mean","Good-Fellows"))
	if(is.na(weight.type)) { # User typed an unrecognizable name
		stop(paste('You must specify a valid weight.type. The valid types are "Gile\'s SS","RDS-I", "RDS-II", "RDS-I (DS)", "Good-Fellows" and "Arithmetic Mean"'), call.=FALSE)
	}
	if(weight.type %in% c("Gile's SS","Good-Fellows") && is.null(N)){
		stop(paste(weight.type,"estimator requires an estimated population size (N)"))
	}
	
	weights <- switch(weight.type, 
			`RDS-I` = rds.I.weights(rds.data,N=N,...), 
			`RDS-I (DS)` = rds.I.weights(rds.data,N=N,smoothed=TRUE,...),
			`RDS-II` = vh.weights(degs = deg, N=N),
			`Arithmetic Mean` = rep(ifelse(is.null(N),1,N)/n, n),
			`Gile's SS` = gile.ss.weights(degs = deg, N = N, SS.infinity=control$SS.infinity, ...)
	)
#			`Good-Fellows` = gf.weights(rds.data, N = N, SS.infinity=control$SS.infinity, ...)
	se <- substitute(subset)
	if(!is.null(se)){
		if(class(se)!="name") subset <- eval(subset, rds.data, parent.frame())
		subset[is.na(subset)] <- FALSE
		a <- weights[subset]
		if(weight.type=="Gile's SS"){
			attr(a,"N") <- attr(weights,"N")
			attr(a,"estimateN") <- attr(weights,"estimateN")
		}
		weights <- a
	}
	return(weights)
}



#' RDS-I weights
#' @param rds.data An rds.data.frame
#' @param outcome.variable The variable used to base the weights on.
#' @param N Population size
#' @param smoothed Should the data smoothed RDS-I weights be computed.
#' @param ... Unused
#' @export
rds.I.weights<-function(rds.data, outcome.variable, N=NULL,smoothed=FALSE,...){
	if(is.null(rds.data[[outcome.variable]])){
		stop(paste("RDS-I outcome.variable", outcome.variable,"not present in data"))
	}
	tij <- count.transitions(rds.data, outcome.variable)
	network.size <- attr(rds.data, "network.size.variable")
	remvalues <- !is.na(rds.data[[network.size]]) & (rds.data[[network.size]] > 0)
	if(sum(remvalues) < nrow(rds.data)){
		warning(paste(nrow(rds.data)-sum(remvalues),"of",nrow(rds.data),
						"network size values were missing and were removed."), call. = FALSE)
		rds.data[[network.size]][!remvalues] <- max(rds.data[[network.size]],na.rm=TRUE)+1
	}
	if(ncol(tij)>1){
		smoothed.tij <- matrix(nrow=0,ncol=0)
		markov.mle <- prop.table(tij, margin=1)
		q.hat <- get.stationary.distribution(markov.mle)
		if(smoothed){
#			Demographic adjustment of raw transition counts
			smoothed.tij <- sum(tij)*q.hat*markov.mle
			smoothed.tij <- 0.5 * (smoothed.tij + t(smoothed.tij))
			markov.mle <- prop.table(smoothed.tij, margin=1)
			q.hat <- get.stationary.distribution(markov.mle)
		}
		h.hat <- get.h.hat(rds.data,outcome.variable,network.size)    
		est <- as.list((q.hat/h.hat)/sum(q.hat/h.hat))
		
		d=tapply(rds.data[[outcome.variable]],rds.data[[outcome.variable]],length)
		f=match(rds.data[[outcome.variable]],names(est))
		weights <- rep(0,length(f))
		weights[!is.na(f)]=as.numeric(unlist(est[f]))
		weights <- as.numeric(weights/d[f])
		weights[is.na(weights)] <- 0
	}else{
		weights <- rep(1/length(rds.data[[outcome.variable]]),length(rds.data[[outcome.variable]]))
	}
	if(!is.null(N)){
		weights <- N * weights / sum(weights)
	}
	weights
}


#' Volz-Heckathorn (RDS-II) weights
#' @param degs The degrees (i.e. network sizes) of the sample units.
#' @param N Population size
#' @export
vh.weights<-function(degs, N=NULL){
	if(is.null(degs)){
		return(NULL)
	}
	degs[degs==0]<-1
	isnadegs <- is.na(degs)
	degs <- sum(!isnadegs)*(degs)/sum(degs,na.rm=TRUE)
	weights <- (1/degs)
	weights[is.na(weights)] <- 0
	if(!is.null(N)){
		weights <- N * weights / sum(weights)
	}
	weights
}


#was spps_weights
#' Weights using Giles SS estimator
#' @param degs subjects' degrees (i.e. network sizes).
#' @param N Population size estimate.
#' @param number.ss.samples.per.iteration The number of samples to use to estimate inclusion probabilities
#' in a probability proportional to size without replacement design.
#' @param number.ss.iterations number of iterations to use in giles SS algorithm.
#' @param hajek Should the hajek estiamtor be used. If false, the HT estimator is used.
#' @param SS.infinity The sample proportion, \code{n/N}, below which the computation of the SS weights should simplify to that of the \code{RDS-II} weights.
#' @param se Should covariances be included.
#' @param ... unused
#' @export
gile.ss.weights<-function(degs,N,number.ss.samples.per.iteration=500,number.ss.iterations=5,
		hajek=TRUE,SS.infinity=0.04,se=FALSE,...){
	if(is.null(degs)){
		return(NULL)
	}
	degs[degs==0]<-1
	isnadegs <- is.na(degs)
	n <- length(degs[!isnadegs])
	if(n==0){
		return(NULL)
	}
	if(n==1){
		weights <- rep(0,length(degs))
		weights[!isnadegs] <- N
		return(weights)
	}
	mapping<-getestCstacked(samp=degs[!isnadegs],n=N,nit=number.ss.iterations,
			nsampsamp=number.ss.samples.per.iteration,trace=FALSE,hajek=hajek,SS.infinity=SS.infinity)
	if(!hajek){
		sprintf("The estimated population size is %d.\n",mapping$n)
	}
	weights <- degs
	pis=stats::approx(x=mapping$classes,y=1/mapping$probs,xout=degs[!isnadegs],rule=2)$y
	weights[!isnadegs] <- pis
	weights[is.na(weights)] <- 0
	weights <- mapping$n*weights/sum(weights,na.rm=TRUE)
	attr(weights,"N") <- N
	attr(weights,"estimateN") <- mapping$n
	if(se){
		# Set up for covariances
		theta <- mapping$props
		covtheta <- -outer(theta, theta) / mapping$n
		diag(covtheta) <- (1-theta)*theta / mapping$n
		cov <- covtheta
		cov[is.na(cov)] <- 0
		cov <- cov / outer(mapping$n*theta*mapping$probs, mapping$n*theta*mapping$probs)
		#   covariances of the (normalized) weights
		md <- match(degs,mapping$classes)
		wcov <- outer(md,md,function(i,j){cov[cbind(i,j)]})
		attr(weights,"cov") <- wcov
		attr(weights,"classes") <- degs
	}
	weights
}
