# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date: December 2009
# Version 0.1
# Licence GPL v3

setClass('DistModel',
	contains = 'VIRTUAL',
	representation (
		presence = 'data.frame',
		absence = 'data.frame',
		hasabsence = 'logical'
	),	
	prototype (	
		presence = data.frame(),
		absence = data.frame(),
		hasabsence = FALSE
	),
	validity = function(object)	{
		if (object@hasabsence) {
			t1 <- ncol(object@presence) == ncol(object@absence)
			t2 <- sort(colnames(object@presence)) == sort(colnames(object@absence))
			return(t1 & t2)
		} else {
			return(TRUE)
		} 
	}
)	



setMethod ('show' , 'DistModel', 
	function(object) {
		cat('class    :' , class(object), '\n\n')
		cat('variables:', colnames(object@presence), '\n\n')
		pp <- nrow(object@presence)
		cat('\npresence points:', pp, '\n')
		if (pp < 10) {
			print(object@presence)
		} else {
			print(object@presence[1:10,])
			cat('  (... ...  ...)\n\n')
		}
		if (object@hasabsence) {
			pp <- nrow(object@absence)
			cat('\nabsence points:', pp, '\n')
			if (pp < 10) {
				print(object@absence)
			} else {
				print(object@absence[1:10,])
				cat('  (... ...  ...)\n\n')
			}
		}
	}
)	



setClass('ModelEvaluation',
	representation (
		presence = 'vector',
		absence = 'vector',
		np = 'integer',
		na = 'integer',
		auc = 'numeric',
		pauc = 'numeric',
		cor = 'numeric',
		pcor = 'numeric',
		t = 'vector',
		confusion = 'matrix',
		prevalence = 'vector',
		ODP = 'vector', # overall diagnostic power
		CCR = 'vector', # correct classification rate
		TPR = 'vector', # sensitivity, or true positive rate
		TNR = 'vector', # specificity, or true negative rate
		FPR ='vector',  # False positive rate
		FNR ='vector',  # False negative rate
		PPP = 'vector',
		NPP = 'vector',
		MCR = 'vector', # misclassification rate
		OR = 'vector',  # odds ratio
		kappa = 'vector'
	),	
	prototype (	
		np = as.integer(0),
		na = as.integer(0)
	),
	validity = function(object)	{
		return(TRUE)
	}
)

