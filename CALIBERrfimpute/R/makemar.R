makemar <-
function(simdata, prop=0.2){
	# generates a probability of missingness for x1 and x2 which is
	# based on the logistic of y + x3 (i.e. it is dependent on outcome and
	# fully observed)

	predictions <- function(lp, n){
		# uses the vector of linear predictions (lp) from a logistic model
		# and the expected number of positive responses (n) to generate
		# a set of predictions by modifying the baseline
		# n does not have to be an integer
		logistic <- function(x){
			exp(x)/(1+exp(x))
		}
		trialn <- function(lptrial){
			sum(logistic(lptrial))
		}
		stepsize <- 32
		lptrial <- lp
		while(abs(trialn(lptrial) - n) > 1){
			if (trialn(lptrial) > n){
				# trialn bigger than required 
				lptrial <- lptrial - stepsize
			} else {
				lptrial <- lptrial + stepsize
			}
			stepsize <- stepsize / 2
		}
		# Generate predictions from binomial distribution
		pred <- as.logical(rbinom(logical(length(lp)), 1, logistic(lptrial)))
		list(offset=(lptrial-lp)[1], pred=pred, n=sum(pred))
	}

	simdata[predictions(simdata[,'y'] + simdata[,'x3'],
		prop*nrow(simdata))$pred, 'x1'] <- NA
	simdata[predictions(simdata[,'y'] + simdata[,'x3'],
		prop*nrow(simdata))$pred, 'x2'] <- NA
	return(simdata)	
}
