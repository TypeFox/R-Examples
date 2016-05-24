
differential.activity.estimates.local <- function(rds.data, 
		outcome.variable, weight.type = "Gile's SS", uncertainty = "Gile's SS", 
		N, subset, ...) {
	
	if(is(rds.data,"rds.data.frame")){
		
		if(!(outcome.variable %in% names(rds.data))){
			stop(sprintf("No variable called %s appears in the data.",
							outcome.variable))
		}
		network.size <- attr(rds.data, "network.size.variable")
		
	}
	else{
		stop("rds.data must be of type rds.data.frame")
	}
	
	if(weight.type %in% c("RDS-I","RDS-I (DS)")){
		outcome <- factor(rds.data[[outcome.variable]])
# 		Make sure the factor labels are alphabetic!
		outcome=factor(outcome,levels=levels(outcome)[order(levels(outcome))])
		rds.data[[outcome.variable]] <- outcome
	}
	
	
	if(is.null(N)){
		N <- attr(rds.data, "population.size.mid")
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
	
	#########################################################################################
	# The simple data management tasks have been taken care of, so it's now time to compute 
	# the RDS estimates.   The cases for numeric and categorical outcomes are handled 
	# separately.
	
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
	weights.nomiss <- compute.weights(rds.data.nomiss,
			weight.type=weight.type,outcome.variable=outcome.variable, N=N, ...)
	
	outcome <- as.vector(rds.data.nomiss[[outcome.variable]])
	toutcome <- table(outcome)
	maxa <- which.max(toutcome)
	toutcome[maxa] <- 0
	maxa <- c(maxa, which.max(toutcome))
	onames <- sort(names(toutcome)[maxa])
	outcome <- as.numeric(outcome == max(onames, na.rm = TRUE))
	#outcome[is.na(outcome)] <- 0
	moutcome <- max(outcome, na.rm = TRUE)
	prop.outcome <- HT.estimate(weights = weights.nomiss, outcome = outcome)
	deg <- get.net.size(rds.data.nomiss)
	mean.degree.outcome0 <- HT.estimate(weights = weights.nomiss, 
			outcome = deg * (outcome == 0))/(1 - prop.outcome)
	mean.degree.outcome1 <- HT.estimate(weights = weights.nomiss, 
			outcome = deg * (outcome == 1))/prop.outcome
	
	result <- c(mean.degree.outcome1/mean.degree.outcome0, moutcome)
	class(result) <- "differential.activity.estimate"
	result
}

#' Differential Activity between groups
#' @param rds.data An rds.data.frame object
#' @param outcome.variable A character string of column names representing categorical variables.
#' @param weight.type A string giving the type of estimator to use. The options
#' are \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I/DS"},
#' and \code{"Arithemic Mean"}. It defaults to \code{"Gile's
#' SS"}.
#' @param N The population size.
#' @param subset An expression defining a subset of rds.data.
#' @param ... Additional parameters passed to compute.weights.
#' @details This function estimates the ratio of the average degree of one population
#' group divided by the average degree of those in another population group.
#' @examples 
#' data(faux)
#' differential.activity.estimates(faux,"X",weight.type="RDS-II")
#' @export
differential.activity.estimates <- function(rds.data, 
		outcome.variable, weight.type = "Gile's SS", 
		N = NULL, subset = NULL, ...) {
	se <- substitute(subset)
	if(is.null(se)){
	  csubset <- ""
	}else{
	  csubset <- as.character(enquote(substitute(subset)))[2]
	}
	subset <- eval(se,rds.data,parent.frame())
	if (length(outcome.variable) == 1) {
		result <- differential.activity.estimates.local(rds.data=rds.data, 
						outcome.variable=outcome.variable, weight.type = weight.type, 
						N = N, subset = subset, ...)
	}
	else {
		result <- lapply(X = outcome.variable, FUN = function(g) {
					differential.activity.estimates.local(rds.data=rds.data, 
							outcome.variable=g, weight.type = weight.type, 
							N = N, subset = subset, ...)
				})
		names(result) <- outcome.variable
		
	}
	result
}


#' Prints an differential.activity.estimate object
#' @param x an differential.activity.estimate object
#' @param ... unused
#' @method print differential.activity.estimate
#' @export
print.differential.activity.estimate <- function(x,...){
		cat(paste("The mean degree of those with value", x[2], 
						"divided by the mean degree of those without is"), 
				format(x[1], width = 6), "\n")
}
