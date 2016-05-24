summary.WCE <- function(object, allres = FALSE, ...){
	objname <- deparse(substitute(object))
	if(allres == TRUE) {
		.sumWCEall(object, objname)
	}
	else if (allres == FALSE) {
		.sumWCEbest(object, objname)
	}
cat('\nIf you report these results, please cite Sylvestre MP, Abrahamowicz M. Flexible Modeling of the Effects of Time-Dependent Exposures on the
Hazard. Statistics in Medicine 2009; 28(27):3437-3453.\n')
}