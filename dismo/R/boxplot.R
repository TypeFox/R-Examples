

if (!isGeneric("boxplot")) {
	setGeneric("boxplot", function(x, ...)
		standardGeneric("boxplot"))
}	


setMethod('boxplot', signature(x='ModelEvaluation'), 
	function(x, notch=TRUE, ...) {
		comb <- c(x@presence, x@absence)
		group <- c(rep('presence', length(x@presence)), rep('absence', length(x@absence)) )
		boxplot(comb~group, notch=notch, data=data.frame(comb,group), ...)
	}
)


