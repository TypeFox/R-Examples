#'Summarizing TcGSA
#'
#'\code{summary} method for class '\code{TcGSA}'
#'
#'
#'@aliases summary.TcGSA print.summary.TcGSA
#'
#'
#'@param object 
#'an object of class '\code{TcGSA}'.
#'
#'@param \dots further arguments passed to or from other methods.
#'
#'@return The function summary.TcGSA returns a list with the following
#'components (list elements):\itemize{
#'\item time_func the chosen form for the time trend.
#'\item separateSubjects a logical flag indicating wether gene sets
#'tested for discriminating among patients, or for time trends over time.
#'\item ntg the number of treatment groups.
#'\item ngs the number of tested gene sets.
#'\item nsignif the number of significant gene sets at a 5\% FDR (using
#'the default Benjamini & Yekutieli step-up procedure).
#'}
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{TcGSA.LR}}
#'
#'@method summary TcGSA
#'
#'@export
#'
#'@examples
#'
#'data(data_simu_TcGSA)
#'
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'summary(tcgsa_sim_1grp)
#'
#'\dontrun{
#'tcgsa_sim_2grp <- TcGSA.LR(expr=expr_2grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE, 
#'                           group_name="group.var")
#'summary(tcgsa_sim_2grp)
#'}
#'
summary.TcGSA <-function(object, ...){
	signifRes <- signifLRT.TcGSA(object, ...)
  nsignif <- dim(signifRes$mixedLRTadjRes)[1]
  time_func <- object[["time_func"]]
  separateSubjects <- object[["separateSubjects"]]
  ntg <- ifelse(is.null(object[["group.var"]]),1,length(levels(object[["group.var"]])))
  ngs <- length(object[["GeneSets_gmt"]]$geneset.name)
  res  <- list("time_func"=time_func, "separateSubjects"=separateSubjects, "ntg"=ntg, "ngs"=ngs, 
  						 "nsignif"=nsignif, "threshold"=signifRes$threshold, "multCorProc"=signifRes$multCorProc)
  class(res) <- "summary.TcGSA"
  
  return(res)
}

#'@rdname summary.TcGSA
#'
#'@param x an object of class '\code{summary.TcGSA}'.
#'
#'@method print summary.TcGSA
#'
#'@export
#'
print.summary.TcGSA <-function(x, ...){
	cat("\t\tA TcGSA object")
	cat("\n")
	cat("Form of the time trend:")
	cat("\n\t")
	cat(x[["time_func"]])
	cat("\n")
	cat("Number of treatment groups:")
	cat("\n\t")
	cat(x[["ntg"]])
	cat("\n")
	if(x[["separateSubjects"]]){
		cat("Number of gene sets tested for discriminating time trends among patients:") 
	}else{
		cat("Number of gene sets tested for significant time trend:")
	}
	cat("\n\t")
	cat(x[["ngs"]])
	cat("\n\n")
	cat("Number of significant gene sets at a ", x$threshold*100,"% threshold (", x$multCorProc, " procedure):", sep="")
	cat("\n\t")
	cat(x[["nsignif"]])
	cat(" out of ")
	cat(x[["ngs"]])
	cat(" gene sets")
}
