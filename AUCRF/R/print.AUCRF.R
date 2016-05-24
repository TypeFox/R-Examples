print.AUCRF <-
function(x,...){
	cat("\nNumber of selected variables: Kopt=",x$Kopt,"\n")
  cat("AUC of selected variables: OOB-AUCopt=",x$"OOB-AUCopt","\n")
  if(!is.null(x$cvAUC)) cat("AUC from cross validation:",x$cvAUC,"\n")
  cat("Importance Measure:",x$ImpMeasure,"\n")
  cat("\n")

	if(is.null(x$Psel))
    res <- data.frame("Selected Variables"=x$Xopt, Importance=x$ranking[x$Xopt])
  else
  	res <- data.frame("Selected Variables"=x$Xopt, Importance=x$ranking[x$Xopt], Prob.Select=x$Psel[x$Xopt])
  rownames(res) <- NULL
  print(res)
}
