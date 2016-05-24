

RDS.bootstrap.intervals.local <- function(rds.data, 
		outcome.variable, weight.type, uncertainty, N, subset,number.of.bootstrap.samples, 
		confidence.level=.95, 
		control=control.rds.estimates(), continuous=NA, fast=FALSE, useC=FALSE, csubset="", ci.type="t", ...) {
	
	if (is(rds.data, "rds.data.frame")) {
		if (!(outcome.variable %in% names(rds.data))) {
			stop(sprintf("No variable called %s appears in the data.", 
							outcome.variable))
		}
		network.size <- attr(rds.data, "network.size.variable")
		if(!(network.size %in% names(rds.data)))
			stop("invalid network size variable")
		else{
	         is.cts <- FALSE
	         is.quantile <- FALSE
#	         if( storage.mode(rds.data[[outcome.variable]])=="double"  |
#	             (!is.na(continuous) && length(continuous)==1 && continuous < 1 && continuous > 0) ){is.cts <- TRUE}
 	         if( !is.na(continuous) && length(continuous)==1 ){
		   if( is.character(continuous) && continuous == "mean" ){is.cts <- TRUE}
 	           if( is.numeric(continuous)   && continuous < 1 && continuous > 0 ){is.cts <- TRUE;is.quantile <- TRUE}
		 }
		 rds.data[[network.size]] <- as.numeric(rds.data[[network.size]])
		}
	}
	else {
		stop("rds.data must be of type rds.data.frame")
	}
	
	if (is.null(weight.type)) {
		weight.type <- "Gile's SS"
	}
	weight.type <- match.arg(weight.type,
			c("Gile's SS","RDS-I", "RDS-II", "RDS-I (DS)","Good-Fellows","Arithmetic Mean"))
	if(is.na(weight.type)) { 
		# User typed an unrecognizable name
		stop(paste('You must specify a valid weight.type. The valid types are "Gile\'s SS","RDS-I", "RDS-II", "RDS-I (DS)", "Good-Fellows", and "Arithmetic Mean"'), call.=FALSE)
	}
	
	if (is.null(uncertainty)) {
		if(weight.type %in% c("Gile's SS","Good-Fellows")){
			uncertainty <- "Gile"
		}else if(weight.type %in% c("RDS-I","RDS-I (DS)","RDS-II")){
			uncertainty <- "Salganik"
		}else if(weight.type=="Arithmetic Mean"){
			uncertainty <- "SRS"
		}
	}
	uncertainty <- match.arg(uncertainty,c("Gile","Salganik","SRS"))
	if (is.null(N)) {
		N <- attr(rds.data, "population.size.mid")
		if (is.null(N) && uncertainty == "Gile"){
			stop("N not specified with no default in the data. N is required for Gile's SS estimator and bootstrap")
		}
		if(!is.null(N))
			cat("\nNote: Using the data's mid population size estimate: N =", 
				N, "\n")
	}
	
	###########################################################################################
	# Check for missing values and warn the user if any are removed.   This should really taken
	# care of elsewhere.  NB: It is also worth considering the semantics of the message 
	# "values were missing and were removed".    
	remvalues <- rds.data[[network.size]]==0 | is.na(rds.data[[network.size]])
	if(any(remvalues)){
		warning(paste(sum(remvalues),"of",nrow(rds.data),
						"network sizes were missing or zero. The estimator will presume these are",max(rds.data[[network.size]],na.rm=TRUE)), call. = FALSE)
		
		rds.data[[network.size]][remvalues] <- max(rds.data[[network.size]],na.rm=TRUE)
	}
	rds.data.nomiss <- rds.data
	
	se <- substitute(subset)
	subset <- eval(se,rds.data,parent.frame())
	if(is.null(se)|is.null(subset)){
		subset <- rep(TRUE,length=nrow(rds.data.nomiss))
	}else{
		subset[is.na(subset)] <- FALSE
		if(!is.null(N)){
			#use VH estimator to adjust population size to sub-population
			tmp.wts <- vh.weights(rds.data.nomiss[[network.size]])
			tmp.wts <- tmp.wts/sum(tmp.wts)
			prop <- sum(tmp.wts*subset)
			N <- N * prop
			if(N < sum(subset))
				stop("Estimated sub-population size smaller than subset.")
		}
		
		#This subsets setting orphaned children to be seeds
		#in order to maintain a valid recruitment tree.
		rds.data.nomiss <- rds.data.nomiss[subset,,warn=FALSE]
		
		#drop 0 count levels
		if(is.factor(rds.data[[outcome.variable]])){
			outcome <- factor(rds.data.nomiss[[outcome.variable]])
# 			Make sure the factor labels are alphabetic!
			outcome=factor(outcome,levels=levels(outcome)[order(levels(outcome))])
			rds.data.nomiss[[outcome.variable]] <- outcome
		}
	}

	if( mean(duplicated(rds.data.nomiss[[outcome.variable]],na.rm=TRUE)) < control$discrete.cutoff ){is.cts <- TRUE}
	
	#only categorical estimates are supported
	outcome <- factor(rds.data.nomiss[[outcome.variable]])
#       Make sure the factor labels are alphabetic!
	outcome <- factor(outcome,levels=levels(outcome)[order(levels(outcome))])
	rds.data.nomiss[[outcome.variable]] <- outcome
	outclasses <- levels(outcome)
	g <- length(outclasses)
	if(g==1){
		warning(paste(outcome.variable,"has only one level. Skipping..."))
		return(invisible(NULL))
	}
	nsamples <- sum(!is.na(as.vector(outcome)))
	nsamplesbyoutcome <- table(as.numeric(outcome))
	
	rds <- RDS.estimates.local(
			rds.data = rds.data.nomiss,
			outcome.variable=outcome.variable,
			N=N,
			weight.type=weight.type,
			control=control,
			continuous=continuous, is.cts=is.cts, is.quantile=is.quantile,
			empir.lik=TRUE, useC=useC,
			...)
	
        if(is.rds.interval.estimate(rds)){
	  observed.estimate <- rds$estimate
	  weights.all <- rds$weights
        }else{
	  observed.estimate <- rds@estimate
	  weights.all <- rds@weights
        }
	
        if(is.null(number.of.bootstrap.samples) || number.of.bootstrap.samples < 10){
          if(is.rds.interval.estimate(rds)){
            EL.se <- matrix(rds$interval, ncol = 6, byrow = FALSE)[1,5]
          }else{
            EL.se <- sqrt(observed.estimate*(1-observed.estimate)/length(weights.all))
          }
          number.of.bootstrap.samples <- min(500,max(20,ceiling(0.5*(EL.se/0.002)^2+1)))
          if(!is.numeric(number.of.bootstrap.samples)){
            number.of.bootstrap.samples <- 500
          }
	  cat(sprintf("Using %d bootstrap samples.\n", number.of.bootstrap.samples))
        }
	#
	#  Confidence intervals
	#
	crit <- qnorm((1-confidence.level)/2,lower.tail=FALSE)
	if (uncertainty == "Salganik") {
		
		bs <- salganik.bootstrap.se(rds.data = rds.data.nomiss, 
				group.variable = outcome.variable, 
				number.of.bootstrap.samples = number.of.bootstrap.samples,
				estimator.name = weight.type, 
				N=N,
				...)

		attr(bs,"is.cts") <- is.cts
		attr(bs,"is.quantile") <- is.quantile
	        attr(bs,"mu") <- observed.estimate
		attr(bs,"sigma2") <- observed.estimate * (1 - observed.estimate)
		
		bso <- bs[match(names(bs),outclasses)]
		estimate <- cbind(observed.estimate, observed.estimate - crit * bso,
				observed.estimate + crit * bso)
		colnames(estimate) <- c("point", "lower", "upper")
		estimate=switch(ci.type,
		"pivotal"={
		  bsests <- attr(bs,"bsresult")
		  bsests <- bsests[,match(names(bs),outclasses)]
		  qcrit <- t(apply(bsests,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
 		  estimate[,2:3] <- 2*observed.estimate - qcrit
#	  	  estimate[, 1] <- 2*observed.estimate-apply(bsests,2,mean,na.rm=TRUE)
		  estimate[, 1] <- observed.estimate
		  estimate},
		"quantile"={
		  bsests <- attr(bs,"bsresult")
		  bsests <- bsests[,match(names(bs),outclasses)]
		  qcrit <- t(apply(bsests,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
		  estimate[,2:3] <- qcrit[,2:1]
		  estimate[, 1] <- observed.estimate
		  estimate},
		"proportion"={
		  if(all(observed.estimate >= 0) & (any(observed.estimate < 0.15) | any(estimate[,2:3]<0|estimate[,2:3]>1))){
		    bsests <- attr(bs,"bsresult")
		    bsests <- bsests[,match(names(bs),outclasses)]
		    qcrit <- t(apply(bsests,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
		    estimate[,2:3] <- qcrit[,2:1]
		  }
		  estimate[estimate[,2]<0,2] <- 0
		  estimate[estimate[,3]<0,3] <- 0
		  estimate[estimate[,2]>1,2] <- 1
		  estimate[estimate[,3]>1,3] <- 1
		  estimate[, 1] <- observed.estimate
		  estimate},
		{# t (the default)
		  estimate[, 1] <- observed.estimate
		  estimate}
                )
	} else if(uncertainty == "SRS"){
		varsrs <- observed.estimate * (1 - observed.estimate)/nsamples
		if(!is.null(N))
			varsrs <- ((N - nsamples)/(N - 1)) * varsrs
		estimate <- cbind(observed.estimate, observed.estimate - 
						crit * sqrt(varsrs), 
				observed.estimate + crit * sqrt(varsrs))
		colnames(estimate) <- c("point", "lower", "upper")
		bs <- varsrs
		attr(bs,"is.cts") <- is.cts
		attr(bs,"is.quantile") <- is.quantile
	        attr(bs,"mu") <- observed.estimate
	        if(is.cts){
		  oc <- rds.data.nomiss[[outcome.variable]]
	          if(!is.numeric(oc)){oc <- as.numeric(as.character(oc))}
		  fn.mu <- function(x,dis){wtd.mean(dis,x,na.rm=TRUE,normwt=TRUE)}
		  attr(bs,"sigma2") <- fn.mu(weights.all,oc^2)-fn.mu(weights.all,oc)^2
		}else{
		  attr(bs,"sigma2") <- observed.estimate * (1 - observed.estimate)
		}

	} else if(uncertainty == "Gile"){
		bs <- SSBS.estimates(rds.data.nomiss, 
				outcome.variable,  
				confidence.level = confidence.level, 
				N = N,
				weight.type=weight.type, 
				number.of.bootstrap.samples = number.of.bootstrap.samples,
				control=control, continuous=continuous, is.cts=is.cts, is.quantile=is.quantile, fast=fast, useC=useC,
				...)

		if(attr(bs,"is.cts")){
			observed.estimate <- bs[1]
			outclasses <- outcome.variable
                        nsamplesbyoutcome <- nrow(rds.data.nomiss)
		}

		estimate <- bs[match(rownames(bs),outclasses),,drop=FALSE]
		colnames(estimate) <- c("point", "lower", "upper")
		bsests <- attr(bs,"bsresult")$bsests
		bsests <- bsests[,match(rownames(bs),outclasses)]
		estimate=switch(ci.type,
		"pivotal"={
		  qcrit <- t(apply(bsests,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
 		  estimate[,2:3] <- 2*observed.estimate - qcrit
#	  	  estimate[, 1] <- 2*observed.estimate-apply(bsests,2,mean,na.rm=TRUE)
		  estimate[, 1] <- observed.estimate
		  estimate},
		"quantile"={
		  qcrit <- t(apply(bsests,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
		  estimate[,2:3] <- qcrit[,2:1]
		  estimate[, 1] <- observed.estimate
		  estimate},
		"proportion"={
		  if(all(observed.estimate >= 0) & (any(observed.estimate < 0.15) | any(estimate[,2:3]<0|estimate[,2:3]>1))){
		    qcrit <- t(apply(bsests,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
		    estimate[,2:3] <- qcrit[,2:1]
#	  	    estimate[, 1] <- 2*observed.estimate-apply(bsests,2,mean,na.rm=TRUE)
		  }
		  estimate[estimate[,2]<0,2] <- 0
		  estimate[estimate[,3]<0,3] <- 0
		  estimate[estimate[,2]>1,2] <- 1
		  estimate[estimate[,3]>1,3] <- 1
		  estimate[, 1] <- observed.estimate
		  estimate},
		"hmg"={
                  nobs <- min(as.vector(nsamplesbyoutcome))
                  dae <- as.vector(differential.activity.estimates(rds.data = rds.data.nomiss,
                                   outcome.variable=outcome.variable, N=N, weight.type=weight.type))[1]
		  if(nobs == as.vector(nsamplesbyoutcome)[1]){dae <- 1/dae}
                  if(nobs+control$lowprevalence[1]*dae < control$lowprevalence[2]){
#                  Apply the combined Agresti-Coull and the bootstrap-t interval of Mantalos and Zografos (2008)
		   wbse <- (attr(bs,"bsresult")$bsnm*bsests+2)/(attr(bs,"bsresult")$bsnm+4)
		   nm <- sum(weights.all,na.rm=TRUE)^2/sum(weights.all^2,na.rm=TRUE)
	  	   woe <- (nm*observed.estimate+2)/(nm+4)
		   sigmamb <- sqrt(wbse*(1-wbse)/attr(bs,"bsresult")$bsnm)
		   tstarn <- sweep(wbse,2,woe,"-")/sigmamb
	  	   tstarn[is.infinite(tstarn)] <- NA
	  	   qcrit <- t(apply(tstarn,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
	  	   estimate[,2:3] <- woe-qcrit*sqrt(woe*(1-woe)/(nm+4))
		  }
		  #
		  estimate[, 1] <- observed.estimate
		  estimate},
		"acbt"={
#                 Apply the combined Agresti-Coull and the bootstrap-t interval of Mantalos and Zografos (2008)
		  wbse <- (attr(bs,"bsresult")$bsnm*bsests+2)/(attr(bs,"bsresult")$bsnm+4)
		  nm <- sum(weights.all,na.rm=TRUE)^2/sum(weights.all^2,na.rm=TRUE)
	  	  woe <- (nm*observed.estimate+2)/(nm+4)
		  sigmamb <- sqrt(wbse*(1-wbse)/attr(bs,"bsresult")$bsnm)
		  tstarn <- sweep(wbse,2,woe,"-")/sigmamb
	  	  tstarn[is.infinite(tstarn)] <- NA
	  	  qcrit <- t(apply(tstarn,2,stats::quantile,c((confidence.level+1)/2,(1-confidence.level)/2),na.rm=TRUE))
	  	  estimate[,2:3] <- woe-qcrit*sqrt(woe*(1-woe)/(nm+4))
		  #
		  estimate[, 1] <- observed.estimate
		  estimate},
		{# t (normal;the default)
#		  estimate[,2:3] <- estimate[,2:3]-matrix(apply(estimate[,2:3],1,mean),ncol=2,nrow=nrow(estimate)) + observed.estimate
		  estimate[, 1] <- observed.estimate
		  estimate}
                )
	}
	#
	# Design effect (for factor outcomes)
	#
	if(attr(bs,"is.cts")){
	 if(attr(bs,"is.quantile")){
	  varoutcome <- continuous*(1-continuous)/(stats::dnorm(observed.estimate,mean=attr(bs,"mu"),sd=sqrt(attr(bs,"sigma2"))))^2
	 }else{
	  varoutcome <- attr(bs,"sigma2")
	 }
	}else{
	  varoutcome <- estimate[, 1] * (1 - estimate[, 1])
	}
	#       Note the finite sample correction factor
	if(!is.null(N)){
		varsrs <- (((N - nsamples)/(N - 1)) * varoutcome/nsamples)
	}else{
		varsrs <- varoutcome/nsamples
	}

	de <- ((estimate[, 3] - estimate[, 2])/(crit * 2))^2/varsrs
	estimate <- cbind(estimate, de, (estimate[, 3] - estimate[, 2])/(crit * 2), nsamplesbyoutcome)
	colnames(estimate)[c(4, 5, 6)] <- c("Design Effect", "s.e.", "n")
	if(attr(bs,"is.cts")){
	  rownames(estimate) <- outcome.variable
	}else{
	  outclasses <- sort(unique(as.vector(rds.data.nomiss[[outcome.variable]])))
	  outclasses[outclasses=="NA.NA"] <- "NA"
	  rownames(estimate) <- outclasses[1:nrow(estimate)]
	  estimate <- as.numeric(estimate)
	  names(estimate)[1:g] <- outclasses[1:g]
	}
	if(exists("bs"))
		attr(estimate,"bsresult") <- attr(bs,"bsresult")
	
	result <- rds.interval.estimate(estimate, 
			outcome.variable,
			weight.type, 
			uncertainty, 
			weights.all,
			N=N,
			csubset=csubset,
			conf.level=confidence.level)
	attr(result,"bsresult") <- attr(estimate,"bsresult")
	
	return(result)
}






#' RDS Bootstrap Interval Estimates
#' 
#' This function computes an interval estimate for one or more categorical
#' variables. It optionally uses attributes of the RDS data set to determine
#' the type of estimator and type of uncertainty estimate to use.
#' 
# Note that the current estimator only works on binary outcome variables. (I think this cane be removed now--non-binary categorical outcome variables work.)
#' 
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical or numeric variable to be
#' analyzed.
#' @param weight.type A string giving the type of estimator to use. The options
#' are \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I (DS)"},
#' and \code{"Arithemic Mean"}. If \code{NULL} it defaults to \code{"Gile's
#' SS"}.
#' @param uncertainty A string giving the type of uncertainty estimator to use.
#' The options are \code{"SRS"}, \code{"Gile"} and \code{"Salganik"}. This is usually
#' determined by \code{weight.type} to be consistent with the estimator's
#' origins. The estimators RDS-I, RDS-I (DS), and RDS-II default to \code{"Salganik"},  "Arithmetic
#' Mean" defaults to \code{"SRS"} and "Gile's SS" defaults to the \code{"Gile"} bootstrap.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{population.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing it defaults to 1000.
#' @param subset An optional criterion to subset \code{rds.data} by. It is a
#' character string giving an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{"seed > 0"} to
#' exclude seeds. It can be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param confidence.level The confidence level for the confidence intervals.
#' The default is 0.95 for 95\%.
#' @param number.of.bootstrap.samples The number of bootstrap samples to take
#' in estimating the uncertainty of the estimator. If \code{NULL} it defaults
#' to the number necessary to compute the standard error to accuracy 0.001.
#' @param continuous A numerical value between 0 and 1 or the character value \code{"mean"} or \code{NA}. 
#' If \code{"mean"} it will estimate the population mean of the \code{outcome.variable}.
#' If it is a numeric value betwern zero and one, estimate this quantile of the population distribution of the
#' \code{outcome.variable}. Otherwise it will compute the population frequencies of each value of the \code{outcome.variable}.
#' @param fast Use a fast bootstrap where the weights are reused from the
#' estimator rather than being recomputed for each bootstrap sample.
#' @param useC Use a C-level implementation of Gile's bootstrap (rather than
#' the R level). The implementations should be a computational 
#' equivalent estimator (except for speed).
#' @param ci.type Type of confidence interval to use, if possible. If "t", use lower and upper confidence
#' interval values based on the standard deviation of the bootstrapped values and a t multiplier.
#' If "pivotal", use lower and upper confidence interval values based on the 
#' basic bootstrap (also called the pivotal confidence interval).
#' If "quantile", use lower and upper confidence interval values based on the 
#' quantiles of the bootstrap sample.
#' If "proportion", use the "t" unless the estimated proportion is less than 0.15 or the bounds are 
#' outside [0,1 . In this case, try the "quantile" and constrain the bounds to be compatible with [0,1]. 
#' @param control A list of control parameters for algorithm
#' tuning. Constructed using\cr
#' \code{\link{control.rds.estimates}}.
#' @param ... Additional arguments for RDS.*.estimates.
#' @return An object of class \code{rds.interval.estimate} summarizing the inference.
#' The confidence interval and standard error are based on the bootstrap procedure.
#' In additon, the object has attribute \code{bsresult} which provides details of the 
#' bootstrap procedure. The contents of the \code{bsresult} attribute depends on the
#' \code{uncertainty} used. If \code{uncertainty=="Salganik"} then \code{bsresult} is a
#' vector of standard deviations of the bootstrap samples. 
#' If \code{uncertainty=="Gile's SS"} then
#' \code{bsresult} is a list with components for the bootstrap point estimate,
#' the bootstrap
#' samples themselves and the standard deviations of the bootstrap samples. 
#' If \code{uncertainty=="SRS"} then \code{bsresult} is NULL. 
#' @references Gile, Krista J. 2011 \emph{Improved Inference for
#' Respondent-Driven Sampling Data with Application to HIV Prevalence
#' Estimation}, \emph{Journal of the American Statistical Association}, 106,
#' 135-146.
#' 
#' Gile, Krista J., Handcock, Mark S., 2010 \emph{Respondent-driven Sampling:
#' An Assessment of Current Methodology}. Sociological Methodology 40, 285-327.
#' 
#' @keywords survey manip
#' @examples
#' 
#' \dontrun{
#' data(fauxmadrona)
#' RDS.bootstrap.intervals(rds.data=fauxmadrona,weight.type="RDS-II",
#'      uncertainty="Salganik",
#' 	outcome.variable="disease",N=1000,number.of.bootstrap.samples=50)
#' }
#' 
#' @export
RDS.bootstrap.intervals <- function(rds.data, outcome.variable, 
		weight.type = NULL, uncertainty = NULL, N = NULL, subset = NULL, 
		confidence.level = .95, number.of.bootstrap.samples = NULL, continuous=NA, fast=TRUE, useC=TRUE,
		ci.type="t", control=control.rds.estimates(),
		...) {
	se <- substitute(subset)
	if(is.null(se)){
	  csubset <- ""
	}else{
	  csubset <- as.character(enquote(substitute(subset)))[2]
	}
	subset <- eval(se,rds.data,parent.frame())
	if (length(outcome.variable) == 1) {
		result <- RDS.bootstrap.intervals.local(rds.data, outcome.variable, 
				weight.type, uncertainty, N, subset, confidence.level, 
				number.of.bootstrap.samples=number.of.bootstrap.samples,
				control=control,continuous=continuous,fast=fast,useC=useC,ci.type=ci.type,csubset=csubset )
	}
	else {
		result <- lapply(X = outcome.variable, FUN = function(g) {
					RDS.bootstrap.intervals.local(rds.data, g, weight.type, 
							uncertainty, N, subset, confidence.level,
							number.of.bootstrap.samples=number.of.bootstrap.samples, 
							control=control,continuous=continuous,fast=fast, useC=useC, ci.type=ci.type,							   csubset=csubset, ...)
				})
		names(result) <- outcome.variable
	}
	return(result)
}
