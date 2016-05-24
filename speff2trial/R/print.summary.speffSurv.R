print.summary.speffSurv <- function(x,...){
	if (!x$fixed){
		cat("\nOptimal models using",x$method,"method:\n")
		cat("Rnd space:   ",format(x$formula$rndSpace[-2]),"\n",sep="")
		cat("Cens space:  ",format(x$formula$censSpace[-2]),"\n",sep="")
	}
	cat("\nTreatment effect\n")
	print(x$tab, digits=5, print.gap=2)
}
