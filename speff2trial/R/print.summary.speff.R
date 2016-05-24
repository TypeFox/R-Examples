print.summary.speff <- function(x,...){
	if (!is.null(x$rsq)){
		cat("\nOptimal models using",x$method,"method:\n")
		if (x$predicted[1]) cat("Control:   ",format(x$formula$control),", R-squared: ",round(x$rsq[1],2),"\n",sep="")
		if (x$predicted[2]) cat("Treatment:  ",format(x$formula$treatment),", R-squared:  ",round(x$rsq[2],2),"\n",sep="")
	}
	cat("\nTreatment effect\n")
	print(x$tab, digits=5, print.gap=2)
}
