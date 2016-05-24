#     ____________________________
#    |\_________________________/|\
#    ||                         || \
#    ||    algo.farrington      ||  \
#    ||     new version         ||  |
#    ||                         ||  |
#    ||                         ||  |
#    ||                         ||  |
#    ||                         ||  |
#    ||                         ||  /
#    ||_________________________|| /
#    |/_________________________\|/
#       __\_________________/__/|_
#      |_______________________|/ )
#    ________________________    (__
#   /oooo  oooo  oooo  oooo /|   _  )_
#  /ooooooooooooooooooooooo/ /  (_)_(_)
# /ooooooooooooooooooooooo/ /    (o o)
#/C=_____________________/_/    ==\o/==

# Version of the 26.06.2013
# M.Salmon, M.Hoehle

################################################################################
# CONTENTS
################################################################################
# # MAIN FUNCTION
# Function that manages input and output.
# # RESIDUALS FUNCTION
# Function that calculates Anscombe residuals.
# # WEIGHTS FUNCTION
# Function that calculates weights based on these residuals.
# # FORMULA FUNCTION
# Function that writes a formula for the glm using Booleans from control.
# # FIT GLM FUNCTION
# Function that fits a GLM. If it does not converge this function tries to fit it without time trend.
# # THRESHOLD FUNCTION
# Function that calculates the lower and upper threshold, the probability of observing a count that is >= observed, and the score.
# There are two versions of this function depending on the method chosen.
# # BLOCKS FUNCTION
# Function that creates the factor variable for the glm.
# # DATA GLM FUNCTION
# Function that prepares data for the glm
# # GLM FUNCTION
# Function that calls fit glm, checkst he time trend and calculate the prediction fort he current timepoint.

################################################################################
# END OF CONTENTS
################################################################################


################################################################################
# MAIN FUNCTION
################################################################################

farringtonFlexible <- function(sts, control = list(range = NULL, b = 3, w = 3,
													reweight = TRUE,
													weightsThreshold = 2.58,
													verbose = FALSE,glmWarnings = TRUE,
													alpha = 0.01, trend = TRUE,
													pThresholdTrend = 0.05, limit54=c(5,4),
													powertrans="2/3",
													fitFun="algo.farrington.fitGLM.flexible",
													populationOffset = FALSE,
													noPeriods = 1, pastWeeksNotIncluded = 26,
													thresholdMethod = "delta")) {

    ######################################################################
    # Use special Date class mechanism to find reference months/weeks/days
    ######################################################################

    if (is.null( sts@epochAsDate)) {
        epochAsDate <- FALSE
    } else {
        epochAsDate <-    sts@epochAsDate
    }

    ######################################################################
    # Fetch observed and population
    ######################################################################

    # Fetch observed
    observed <- observed(sts)
    freq <- sts@freq
    if (epochAsDate) {
        epochStr <- switch( as.character(freq), "12" = "month","52" =    "week",
                                                        "365" = "day")
    } else {
        epochStr <- "none"
    }

    # Fetch population (if it exists)
    if (!is.null(population(sts))) {
        population <- population(sts)
    } else {
        population <- rep(1,length(observed))
    }

    ######################################################################
    # Fix missing control options
    ######################################################################
    # How many years to go back in time?
    if (is.null(control[["b",exact=TRUE]])) { control$b = 5 }

    # Half-window length
    if (is.null(control[["w", exact = TRUE]])) { control$w = 3 }

    # Range of time points to be evaluated
    if (is.null(control[["range", exact=TRUE]])) {
        control$range <- (freq*(control$b)+control$w +1):dim(observed)[1]
    }



    # Reweighting past outbreaks?
    if (is.null(control[["reweight",exact=TRUE]])) {control$reweight=TRUE}
    # With which threshold?
    if (is.null(control[["weightsThreshold",exact=TRUE]])) {
        control$weightsThreshold=2.58
    }

    # Printing information?
    if (is.null(control[["verbose",exact=TRUE]]))    {control$verbose=FALSE}

    # Printing warning from glm.fit?
    if (is.null(control[["glmWarnings",exact=TRUE]]))    {control$glmWarnings=TRUE}

    # An approximate (two-sided) (1 - alpha)% prediction interval is calculated
    if (is.null(control[["alpha",exact=TRUE]]))        {control$alpha=0.05}

    # Include a time trend when possible?
    if (is.null(control[["trend",exact=TRUE]]))        {control$trend=TRUE}

    # Which pvalue for the time trend to be significant?
    if (is.null(control[["pThresholdTrend",exact=TRUE]])){
        control$pThresholdTrend=0.05}


    # No alarm is sounded
    #    if fewer than cases = 5 reports were received in the past period = 4
    #    weeks. limit54=c(cases,period) is a vector allowing the user to change
    #    these numbers
    if (is.null(control[["limit54",exact=TRUE]]))    {control$limit54=c(5,4)}

    # Power transformation to apply to the data.
    if (is.null(control[["powertrans",exact=TRUE]])){control$powertrans="2/3"}

    # How many noPeriods between windows around reference weeks?
    if (is.null(control[["noPeriods",exact=TRUE]])){control$noPeriods=1}

    # Use factors in the model? Depends on noPeriods, no input from the user.
    if (control$noPeriods!=1) {
        control$factorsBool=TRUE
    } else {
        control$factorsBool=FALSE
    }

    # Use a population offset in the model?
    if (is.null(control[["populationOffset",exact=TRUE]])) {
        control$populationOffset=FALSE
    }

    # How many past weeks not to take into account?
    if (is.null(control[["pastWeeksNotIncluded",exact=TRUE]])) {
        control$pastWeeksNotIncluded=control$w
    }

    # Which function to use?
    # Only one possibility at the moment
    if (is.null(control[["fitFun",exact=TRUE]])) {
        control$fitFun="algo.farrington.fitGLM.flexible"
    } else {
        control$fitFun <- match.arg(control$fitFun, c(
            "algo.farrington.fitGLM.flexible"))
    }

    # Which method for calculating the threshold?
	# Extracting the method

    if (is.null(control[["thresholdMethod",exact=TRUE]])) { control$thresholdMethod="delta"}
    thresholdMethod<- match.arg(control$thresholdMethod, c("delta","nbPlugin","muan"),several.ok=FALSE)

    # Adapt the argument for the glm function
    control$typePred <- switch(thresholdMethod, "delta"="response","nbPlugin"="link","muan"="link")


    # Which threshold function?
    control$thresholdFunction <- switch(thresholdMethod,
	                                     "delta"="algo.farrington.threshold.farrington",
										 "nbPlugin"="algo.farrington.threshold.noufaily",
										 "muan"="algo.farrington.threshold.noufaily")

    ######################################################################
    # Check options
    ######################################################################
    if (!((control$limit54[1] >= 0) && (control$limit54[2] > 0))) {
        stop("The limit54 arguments are out of bounds: cases >= 0 and period > 0.")
    }

    ######################################################################
    # Initialize the necessary vectors
    ######################################################################
    # Vector for score
    score <- matrix(data = 0, nrow = length(control$range), ncol = ncol(sts))
    sts@control$score <- score

    # Vector for time trend
    trend <- matrix(data = 0, nrow = length(control$range), ncol = ncol(sts))


    # Vector for predictive distribution
    pvalue <- matrix(data = 0, nrow = length(control$range), ncol = ncol(sts))
    sts@control$pvalue <- pvalue

    # Vector for expected count
    expected <- matrix(data = 0, nrow = length(control$range), ncol = ncol(sts))
    sts@control$expected <- expected

    # Vector for mu0 (prediction from glm)
    mu0Vector <- matrix(data = 0, nrow = length(control$range), ncol = ncol(sts))
    sts@control$mu0Vector <- mu0Vector

    # Vector for overdispersion phi (from glm)
    phiVector <- matrix(data = 0, nrow = length(control$range), ncol = ncol(sts))
    sts@control$phiVector <- phiVector

    # Vector for time trend (from glm)
    trendVector <- matrix(data = 0, nrow = length(control$range), ncol = ncol(sts))
    sts@control$trendVector <- trendVector

    # Define objects
    n <- control$b*(2*control$w+1)


    # loop over columns of sts
    for (j in 1:ncol(sts)) {

    	#Vector of dates
		if (epochAsDate) {
			vectorOfDates <- as.Date(sts@epoch, origin="1970-01-01")
		} else {
			vectorOfDates <- seq_len(length(observed[,j]))
		}

        # Loop over control$range
        for (k in control$range) {

            ######################################################################
            # Prepare data for the glm
            ######################################################################
			dayToConsider <- vectorOfDates[k]
			diffDates <- diff(vectorOfDates)
			dataGLM <- algo.farrington.data.glm(dayToConsider=dayToConsider,
			                         b=control$b, freq=freq,
                                     epochAsDate=epochAsDate,
									 epochStr=epochStr,
									 vectorOfDates=vectorOfDates,w=control$w,
									 noPeriods=control$noPeriods,
									 observed=observed[,j],population=population,
									 verbose=control$verbose,
									 pastWeeksNotIncluded=control$pastWeeksNotIncluded,k)

            ######################################################################
            # Fit the model
            ######################################################################
			finalModel <- algo.farrington.glm(dataGLM,timeTrend=control$trend,populationOffset=control$populationOffset,
			                    factorsBool=control$factorsBool,reweight=control$reweight,
                                weightsThreshold=control$weightsThreshold,
								pThresholdTrend=control$pThresholdTrend,b=control$b,
								noPeriods=control$noPeriods,typePred=control$typePred,
								fitFun=control$fitFun,glmWarnings=control$glmWarnings,
								epochAsDate=epochAsDate,dayToConsider=dayToConsider,
								diffDates=diffDates,populationNow=population[k,j],k,
								verbose=control$verbose)
            if (is.null(finalModel)) {
				#Do we have an alarm -- i.e. is observation beyond CI??
				#upperbound only relevant if we can have an alarm (enoughCases)
				sts@alarm[k,j] <- NA
				sts@upperbound[k,j] <- NA

				mu0Vector[(k-min(control$range)+1),j] <- NA
				# Get overdispersion
				phiVector[(k-min(control$range)+1),j] <- NA

				# Get score
				score[(k-min(control$range)+1),j] <- NA

				#Compute bounds of the predictive
				pvalue[(k-min(control$range)+1),j] <- NA

				# Time trend
				trendVector[(k-min(control$range)+1),j] <- NA
				trend[(k-min(control$range)+1),j] <- NA
				warning(paste("The model could not converge with nor without time trend at timepoint ", k," so no result can be given for timepoint ", k,".\n"))
			} else {
				pred <- finalModel$pred
				doTrend <- finalModel$doTrend
				coeffTime <- finalModel$coeffTime




				######################################################################
				# Calculate lower and upper threshold
				######################################################################
				argumentsThreshold <- list(predFit=pred$fit,predSeFit=pred$se.fit,
														phi=finalModel$phi,
															skewness.transform=control$powertrans,
															alpha=control$alpha, y=observed[k,j],
								method=control$thresholdMethod
															)

				lu <- do.call(control$thresholdFunction, args=argumentsThreshold)

				######################################################################
				# Postprocessing steps & output
				######################################################################

				#Compute exceedance score unless less than 5 reports during last 4 weeks.
				#Changed in version 0.9-7 - current week is included now
				enoughCases <- (sum(observed[(k-control$limit54[2]+1):k,j])
												>=control$limit54[1])

				#18 May 2006: Bug/unexpected feature found by Y. Le Strat.
				#the okHistory variable meant to protect against zero count problems,
				#but instead it resulted in exceedance score == 0 for low counts.
				#Now removed to be concordant with the Farrington 1996 paper.
				X <- ifelse(enoughCases,lu$score,NA)

				#Do we have an alarm -- i.e. is observation beyond CI??
				#upperbound only relevant if we can have an alarm (enoughCases)
				sts@alarm[k,j] <- !is.na(X) && (X>1) && observed[k,j]!=0
				sts@upperbound[k,j] <- ifelse(enoughCases,lu$upper,NA)

				# Possible bug alarm although upperbound <- 0?

				# Calculate expected value from glm
				if (is.na(lu$upper)==FALSE) {
					if ( control$typePred=="response"){
						expected[(k-min(control$range)+1),j] <- ifelse(enoughCases,pred$fit,NA)
					} else{
						expected[(k-min(control$range)+1),j] <- ifelse(enoughCases,exp(pred$fit),NA)
					}
				} else {
					expected[(k-min(control$range)+1),j] <- NA
				}

				# Calculate mean of the negbin distribution of the observation
				# Use linear predictor mean and sd
				eta0 <- pred$fit
				seEta0 <- pred$se.fit

				# deduce the quantile for mu0 from eta0 which is normally distributed
				if (control$thresholdMethod=='nbPlugin'){
					mu0Vector[(k-min(control$range)+1),j] <- exp(eta0)
				} else {
					mu0Vector[(k-min(control$range)+1),j] <- exp(qnorm(1-control$alpha, mean=eta0, sd=seEta0))
				}

				# Get overdispersion
				phiVector[(k-min(control$range)+1),j] <- finalModel$phi

				# Get score
				score[(k-min(control$range)+1),j] <- lu$score

				#Compute bounds of the predictive
				pvalue[(k-min(control$range)+1),j] <- lu$prob

				# Time trend
				if(doTrend) {
					trendVector[(k-min(control$range)+1),j] <- coeffTime
					trend[(k-min(control$range)+1),j] <- 1

				} else {
					trendVector[(k-min(control$range)+1),j] <- NA
				}
			}
        }#done looping over all time points
    } #end of loop over cols in sts.
    # Add information about score
    sts@control$score[,j]     <- score[,j]
    # Add information about trend
    sts@control$trend[,j]     <- trend[,j]

    #Add information about predictive distribution
    sts@control$pvalue[,j]     <- pvalue[,j]

    # Add information about expected value
    sts@control$expected[,j] <- expected[,j]

    # Add information about mean of the negbin distribution of the observation
    sts@control$mu0Vector[,j] <- mu0Vector[,j]

    # Add information about overdispersion
    sts@control$phiVector[,j] <- phiVector[,j]

    # Add information about time trend
    sts@control$trendVector[,j] <- trendVector[,j]

    #Done
    return(sts[control$range,])
}

################################################################################
# END OF MAIN FUNCTION
################################################################################

################################################################################
# REFERENCE TIME POINTS FUNCTION
################################################################################
algo.farrington.referencetimepoints <- function(dayToConsider,b=control$b,freq=freq,epochAsDate,epochStr){


	if (epochAsDate) {
		referenceTimePoints <- as.Date(seq(as.Date(dayToConsider,
										origin="1970-01-01"),
										length=(b+1), by="-1 year"))
	} else {
		referenceTimePoints <- seq(dayToConsider, length=(b+1),by=-freq)

		if (referenceTimePoints[b+1]<=0){
			warning("Some reference values did not exist (index<1).")
		}

	}

	if (epochStr == "week") {

		# get the date of the Mondays/Tuesdays/etc so that it compares to
		# the reference data
		# (Mondays for Mondays for instance)

		# Vectors of same days near the date (usually the same week)
		# dayToGet
		dayToGet <- as.numeric(format(dayToConsider, "%w"))
		actualDay <- as.numeric(format(referenceTimePoints, "%w"))
		referenceTimePointsA <- referenceTimePoints -
			(actualDay
			 - dayToGet)

		# Find the other "same day", which is either before or after referenceTimePoints

		referenceTimePointsB <- referenceTimePointsA + ifelse(referenceTimePointsA>referenceTimePoints,-7,7)


		# For each year choose the closest Monday/Tuesday/etc
		# The order of referenceTimePoints is NOT important

		AB <- cbind(referenceTimePointsA,referenceTimePointsB)
		ABnumeric <- cbind(as.numeric(referenceTimePointsA),as.numeric(referenceTimePointsB))
		distMatrix <- abs(ABnumeric-as.numeric(referenceTimePoints))
		idx <- (distMatrix[,1]>distMatrix[,2])+1
		referenceTimePoints <- as.Date(AB[cbind(1:dim(AB)[1],idx)],origin="1970-01-01")

	}

	return(referenceTimePoints)
	}
################################################################################
# END OF REFERENCE TIME POINTS FUNCTION
################################################################################


################################################################################
# RESIDUALS FUNCTION
# anscombe.residuals(m,phi)
# is defined in algo_farrington.R
################################################################################


################################################################################
# WEIGHTS FUNCTION
# algo.farrington.assign.weights(s,weightsThreshold)
# is defined in algo_farrington.R
################################################################################


################################################################################
# FORMULA FUNCTION
################################################################################
# Function for writing the good formula depending on timeTrend,
# populationOffset and factorsBool

formulaGLM <- function(populationOffset=FALSE,timeBool=TRUE,factorsBool=FALSE){
    # Description
    # Args:
    #     populationOffset: ---
    # Returns:
    #     Vector of X

    # Smallest formula
    formulaString <- "response ~ 1"

    # With time trend?
    if (timeBool){
    formulaString <- paste(formulaString,"+wtime",sep ="")}

    # With population offset?

    if(populationOffset){
    formulaString <- paste(formulaString,"+offset(log(population))",sep ="")}


    # With factors?
    if(factorsBool){
    formulaString <- paste(formulaString,"+seasgroups",sep ="")}

    # Return formula as a string
    return(formulaString)
}
################################################################################
# END OF FORMULA FUNCTION
################################################################################




################################################################################
# FIT GLM FUNCTION
################################################################################

algo.farrington.fitGLM.flexible <- function(dataGLM,
timeTrend,populationOffset,factorsBool,reweight,weightsThreshold,glmWarnings,verbose,control,...) {

    # Model formula depends on whether to include a time trend or not.

    theModel <- formulaGLM(populationOffset,timeBool=timeTrend,factorsBool)

    # Fit it -- this is slow. An improvement would be to use glm.fit here.
    # This would change the syntax, however.
    if (glmWarnings) {
        model <- glm(formula(theModel),data=dataGLM,family = quasipoisson(link="log"))
    } else {
        model <- suppressWarnings(glm(formula(theModel),data=dataGLM,family = quasipoisson(link="log")))
    }
    #Check convergence - if no convergence we return empty handed.

    if (!model$converged) {
        #Try without time dependence

        if (timeTrend) {
			theModel <- formulaGLM(populationOffset,timeBool=F,factorsBool)
			if (glmWarnings) {
				model <- glm(as.formula(theModel), data=dataGLM,
												family = quasipoisson(link="log"))
			} else {
				model <- suppressWarnings(glm(as.formula(theModel), data=dataGLM,
											family = quasipoisson(link="log")))
			}
			if (verbose) {cat("Warning: No convergence with timeTrend -- trying without.\n")}
        }

        if (!model$converged) {
        if (verbose) {cat("Warning: No convergence in this case.\n")}
        if (verbose) {print(dataGLM[,c("response","wtime"),exact=TRUE])}
        return(NULL)
        }
    }

    #Overdispersion parameter phi

    phi <- max(summary(model)$dispersion,1)

    #In case reweighting using Anscome residuals is requested

    if (reweight) {
        s <- anscombe.residuals(model,phi)
        omega <- algo.farrington.assign.weights(s,weightsThreshold)
        if (glmWarnings) {
        model <- glm(as.formula(theModel),data=dataGLM,
                                        family=quasipoisson(link="log"),
                                        weights=omega)
        } else {
	    model <- suppressWarnings(glm(as.formula(theModel),data=dataGLM,
                                                                        family=quasipoisson(link="log"),
                                                                        weights=omega))
        }

        #Here, the overdispersion often becomes small, so we use the max
        #to ensure we don't operate with quantities less than 1.
        phi <- max(summary(model)$dispersion,1)
    } # end of refit.


    #Add wtime, response and phi to the model
    model$phi <- phi
    model$wtime <- dataGLM$wtime
    model$response <- dataGLM$response
    model$population <- dataGLM$population
    if (reweight) {
	model$weights <- omega
    } else{
    model$weights <- model$weights
    }
    #Done
    return(model)

}
################################################################################
# END OF FIT GLM FUNCTION
################################################################################




################################################################################
# THRESHOLD FUNCTION FARRINGTON
################################################################################
algo.farrington.threshold.farrington <- function(predFit,predSeFit,phi,
                                                  skewness.transform,
												    alpha,y,method){
	#Fetch mu0 and var(mu0) from the prediction object
	mu0 <- predFit

	tau <- phi + (predSeFit^2)/mu0

	#Standard deviation of prediction, i.e. sqrt(var(h(Y_0)-h(\mu_0)))

	switch(skewness.transform,
		"none" = { se <- sqrt(mu0*tau); exponent <- 1},
		"1/2" = { se <- sqrt(1/4*tau); exponent <- 1/2},
		"2/3"    = { se <- sqrt(4/9*mu0^(1/3)*tau); exponent <- 2/3},
		{ stop("No proper exponent in algo.farrington.threshold.")})

	#Note that lu can contain NA's if e.g. (-1.47)^(3/2)

	lu <- sort((mu0^exponent + c(-1,1)*qnorm(1-alpha)*se)^(1/exponent),
		    na.last=FALSE)

	#Ensure that lower bound is non-negative

	lu[1] <- max(0,lu[1],na.rm=TRUE)

	# probability associated to the observed value as quantile
	q <- pnorm( y^(1/exponent) , mean=mu0^exponent, sd=se,lower.tail=FALSE)


	# calculate score
	x <- ifelse(is.na(lu[2])==FALSE,(y - predFit) / (lu[2] - predFit),NA)
	return(list(lower=lu[1],upper=lu[2],prob=q,score=x))

}

################################################################################
# END OF THRESHOLD FUNCTION FARRINGTON
################################################################################



################################################################################
# THRESHOLD FUNCTION NOUFAILY
################################################################################
algo.farrington.threshold.noufaily <- function(predFit,predSeFit,phi,
                                               skewness.transform,
												    alpha,y,method){

	# method of Angela Noufaily with modifications

	# Use linear predictor mean and sd
	eta0 <- predFit
	seEta0 <- predSeFit

	# deduce the quantile for mu0 from eta0 which is normally distributed
	if (method=='nbPlugin'){
		mu0Quantile <- exp(eta0)
	} else {
		mu0Quantile <- exp(qnorm(1-alpha, mean=eta0, sd=seEta0))
	}

	if (mu0Quantile==Inf){
		lu <- c(NA,NA)
		q <- NA
	# else is when the method is "muan"
	} else{
		# Two cases depending on phi value
		if (phi>1){
			lu<-c(qnbinom(alpha/2,mu0Quantile/(phi-1),1/phi),
			qnbinom(1-alpha/2,mu0Quantile/(phi-1),1/phi))
		} else {
			lu<-c(qpois(alpha/2,mu0Quantile),qpois(1-alpha/2,mu0Quantile))
		}
		# cannot be negative
		lu[1]=max(0,lu[1])

		# probability associated to the observed value as quantile
		if (phi!=1){
			q <- pnbinom(q= y-1 ,size=mu0Quantile/(phi-1),prob=1/phi,lower.tail=FALSE)
		} else{
			q <- ppois(y-1,mu0Quantile,lower.tail=FALSE)
		}


	}
	# calculate score
	x <- ifelse(is.na(lu[2])==FALSE,(y - predFit) / (lu[2] - predFit),NA)
	return(list(lower=lu[1],upper=lu[2],prob=q,score=x))

}
################################################################################
# END OF THRESHOLD FUNCTION NOUFAILY
################################################################################




################################################################################
# BLOCKS FUNCTION
################################################################################
blocks <- function(referenceTimePoints,vectorOfDates,freq,dayToConsider,b,w,p,
epochAsDate) {
    ## INPUT
    # freq: are we dealing with daily/weekly/monthly data?

    # b: how many years to go back in time

    # w: half window length around the reference timepoints

    # p: number of noPeriods one wants the year to be split into

    ## VECTOR OF ABSOLUTE NUMBERS
    # Very useful to write the code!

    vectorOfAbsoluteNumbers <- seq_len(length(vectorOfDates))


    # logical vector indicating where the referenceTimePoints
    # are in the vectorOfDates
    referenceTimePointsOrNot <- vectorOfDates %in%    referenceTimePoints

    ## VECTOR OF FACTORS
    vectorOfFactors <- rep(NA,length(vectorOfDates))

    ## SETTING THE FACTORS
    # Current week
    if (epochAsDate==FALSE){
        now <- which(vectorOfDates==dayToConsider)
    } else {
        now <- which(vectorOfDates==as.Date(dayToConsider))
    }



    vectorOfFactors[(now-w):now]    <- p

    # Reference weeks

    referenceWeeks <- rev(as.numeric(
                    vectorOfAbsoluteNumbers[referenceTimePointsOrNot=='TRUE']))

    for (i in 1:b) {

			# reference week

		refWeek <- referenceWeeks[i+1]

		vectorOfFactors[(refWeek-w):(refWeek+w)] <- p

			# The rest is only useful if ones want factors, otherwise only have
			# reference timepoints like in the old algo.farrington

		if (p!=1){
			# Number of time points to be shared between vectors
			period <- referenceWeeks[i] - 2 * w - 1 - refWeek

			# Check that p is not too big
			if (period < (p-(2*w+1))){stop('Number of factors too big!')}

			# Look for the length of blocks

			lengthOfBlocks <- period %/% (p-1)
			rest <- period %% (p-1)

			vectorLengthOfBlocks <- rep(lengthOfBlocks,p-1)

			# share the rest of the Euclidian division among the first blocks

			add <- seq_len(rest)
			vectorLengthOfBlocks[add] <-    vectorLengthOfBlocks[add]+1

			 # slight transformation necessary for the upcoming code with cumsum
			vectorLengthOfBlocks <- c(0,vectorLengthOfBlocks)

			# fill the vector

			for (j in 1:(p-1)) {
				vectorOfFactors[(refWeek+w+1+cumsum(vectorLengthOfBlocks)[j]):
								(refWeek+w+1+cumsum(vectorLengthOfBlocks)[j+1]-1)]<-j
			}
		}
    }

    ## DONE!

    return(vectorOfFactors) #indent
}

################################################################################
# END OF BLOCKS FUNCTION
################################################################################



################################################################################
# DATA GLM FUNCTION
################################################################################
algo.farrington.data.glm <- function(dayToConsider, b, freq,
                                     epochAsDate,epochStr,
									 vectorOfDates,w,noPeriods,
									 observed,population,
									 verbose,pastWeeksNotIncluded,k){


	# Identify reference time points

	# Same date but with one year, two year, etc, lag
	# b+1 because we need to have the current week in the vector
	referenceTimePoints <- algo.farrington.referencetimepoints(dayToConsider,b=b,
														       freq=freq,
															   epochAsDate=epochAsDate,
														       epochStr=epochStr
															   )

	if (sum((vectorOfDates %in% min(referenceTimePoints)) == rep(FALSE,length(vectorOfDates))) == length(vectorOfDates)){
		stop("Some reference values did not exist (index<1).")
		}

	if (verbose) { cat("k=", k,"\n")}

	# Create the blocks for the noPeriods between windows (including windows)
	# If noPeriods=1 this is a way of identifying windows, actually.

	blocks <- blocks(referenceTimePoints,vectorOfDates,epochStr,dayToConsider,
					b,w,noPeriods,epochAsDate)

	# Here add option for not taking the X past weeks into account
	# to avoid adaptation of the model to emerging outbreaks
	blocksID <- blocks
	blocksID[(k-pastWeeksNotIncluded):k] <- NA

	# Extract values for the timepoints of interest only

	blockIndexes <- which(is.na(blocksID)==FALSE)


	# Time

	# if epochAsDate make sure wtime has a 1 increment
	if (epochAsDate){
		wtime <- (as.numeric(vectorOfDates[blockIndexes])-
								as.numeric(vectorOfDates[blockIndexes][1]))/as.numeric(diff(vectorOfDates))[1]
	} else {
		wtime <-     as.numeric(vectorOfDates[blockIndexes])
	}

	# Factors
	seasgroups <- as.factor(blocks[blockIndexes])

	# Observed
	response <- observed[blockIndexes]

	# Population
	pop <- population[blockIndexes]

	if (verbose) { print(response)}

	dataGLM <- data.frame(response=response,wtime=wtime,population=pop,
						seasgroups=seasgroups,vectorOfDates=vectorOfDates[blockIndexes])
	dataGLM <- dataGLM[is.na(dataGLM$response)==FALSE,]
	return(dataGLM)

}

################################################################################
# END OF DATA GLM FUNCTION
################################################################################


################################################################################
# GLM FUNCTION
################################################################################

algo.farrington.glm <- function(dataGLM,timeTrend,populationOffset,factorsBool,
                                reweight,weightsThreshold,pThresholdTrend,b,
								noPeriods,typePred,fitFun,glmWarnings,epochAsDate,
								dayToConsider,diffDates,populationNow,k,verbose) {

	arguments <- list(dataGLM=dataGLM,
					  timeTrend=timeTrend,
					  populationOffset=populationOffset,
					  factorsBool=factorsBool,reweight=reweight,
					  weightsThreshold=weightsThreshold,glmWarnings=glmWarnings,
					  verbose=verbose,control=control)

	model <- do.call(fitFun, args=arguments)

	#Stupid check to pass on NULL values from the algo.farrington.fitGLM proc.
	if (is.null(model)) return(model)

	######################################################################
	#Time trend
	######################################################################

	#Check whether to include time trend, to do this we need to check whether
	#1) wtime is signifcant at the 95lvl
	#2) the predicted value is not larger than any observed value
	#3) the historical data span at least 3 years.
	doTrend <- NULL

	# if model converged with time trend
	if ("wtime" %in% names(coef(model))){

		# get the prediction for k
		if(epochAsDate){
			wtime=(as.numeric(dayToConsider)-as.numeric(dataGLM$vectorOfDates[1]))/as.numeric(diffDates)[1]
			} else {
			wtime <- c(k)
			}
		pred <- predict.glm(model,newdata=data.frame(wtime=wtime,
												 population=populationNow,
												 seasgroups=factor(noPeriods),
												 dispersion=model$phi),se.fit=TRUE,type="response")

		# check if three criterion ok

		 #is the p-value for the trend significant (0.05) level
		significant <- (summary.glm(model)$coefficients["wtime",4] < pThresholdTrend)

		#have to use at least three years of data to allow for a trend
		atLeastThreeYears <- (b>=3)
		#no horrible predictions
		noExtrapolation <- (pred$fit <= max(dataGLM$response,na.rm=T))

		#All 3 criteria have to be met in order to include the trend. Otherwise
		#it is removed. Only necessary to check this if a trend is requested.
		doTrend <- (atLeastThreeYears && significant && noExtrapolation)

		# if not then refit
		if (doTrend==FALSE) {
			arguments$timeTrend=FALSE
			model <- do.call(fitFun, args=arguments)

		}
	} else {

		doTrend <- FALSE
	}


	#done with time trend
	######################################################################

	######################################################################
	# Calculate prediction                                              #
	######################################################################
	#Predict value

	if(epochAsDate){
	wtime=(as.numeric(dayToConsider)-as.numeric(dataGLM$vectorOfDates[1]))/as.numeric(diffDates)[1]
	} else {
	wtime <- c(k)
	}
	pred <- predict.glm(model,newdata=data.frame(wtime=wtime,
											population=populationNow,
											seasgroups=factor(noPeriods),
											dispersion=model$phi),se.fit=TRUE,type=typePred)
	coeffTime=ifelse(doTrend,summary.glm(model)$coefficients["wtime",1],NA)
	finalModel <- list (pred,doTrend,coeffTime,model$phi)
	names(finalModel) <- c("pred","doTrend","coeffTime","phi")
    return(finalModel)

}




################################################################################
# END OF GLM FUNCTION
################################################################################
