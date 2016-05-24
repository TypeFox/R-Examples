# summary class definitions
# Version:       0.1-6
# Date:     2011-02-24
# Author:         F.M.


summary.impprep <- function(object, nNames=10L,...) {
	if (length(object$blocks) > 1) {
		cat("\n",length(object$blocks), " different missing-data patterns identified.\n",sep="")
	}
	for (i in seq(along=ncol(object$blockNames))) {
	  cat("\n", paste("Block ",i ,": ", paste(head(object$blockNames[ ,i],nNames,...),collapse=" "),"\n",sep=""))
		}
	cat("\nCompletely observed variables:",object$compNames,"\n",sep=" ")
}



########################################################################################


summary.imp <- function(object,...) {
	nMis <- sum(object$indMatrix)
	misShare <- round(nMis/length(object$indMatrix)*100,1)
	cat("\n",misShare,"% missing values in the original data\n")
	impMode <- ifelse(object$M==1,"single imputation", "multiple imputation")
	cat("\nImputation variant: ",impMode,"\n")
	if (object$M > 1) {
		cat("\nNumber of multiple imputations: ",object$M,"\n")
	}
	if (!is.null(object$nIter)) {
		cat("\nNumber of iterations between stored imputations: ",object$nIter,"\n")
	}
}