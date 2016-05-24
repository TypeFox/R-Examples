# Author: Robert J. Hijmans
# Date :  October 2012
# Version 1.0
# Licence GPL v3



if (!isGeneric("threshold")) {
	setGeneric("threshold", function(x, ...)
		standardGeneric("threshold"))
}	


setMethod('threshold', signature(x='ModelEvaluation'),
	function(x, stat='', sensitivity=0.9, ...) {
		r <- list()
		# maximum kappa
		r$kappa <- x@t[which.max(x@kappa)]
		# maximum sum of the sensitivity (true positive rate) and specificity (true negative rate)
		r$spec_sens <- x@t[which.max(x@TPR + x@TNR)]
		# no omission
		r$no_omission <- x@t[max(which(x@confusion[, 'fn'] == 0))]
		# etc
		

		# Suggestions by Diego Nieto-Lugilde		
		# equal prevalence
		r$prevalence = x@t[which.min(abs(x@t - x@prevalence[1]))] 
		# equal sensitivity and specificity
		r$equal_sens_spec <- x@t[which.min(abs(x@TPR - x@TNR))]
		# fixed sensitivity
		r$sensitivity <- x@t[which.min(x@TPR > sensitivity)]
		
    # etc		
		
		r <- data.frame(r)
		rownames(r) <- 'thresholds'
		if (stat != '') {
			stopifnot (stat %in% c('', 'kappa', 'spec_sens', 'no_omission',  'prevalence', 'equal_sens_spec', 'sensitivity')) 
			r[, stat]
		} else {
			r
		}
	}
)

