OptimalSet <-
function(object){
	if(is.null(object$Psel))
		out <- data.frame("Name"=object$Xopt,"Importance"=object$ranking[object$Xopt])
  else
  	out <- data.frame("Name"=object$Xopt,"Importance"=object$ranking[object$Xopt], "Prob.Selection"=object$Psel[object$Xopt])
  
  rownames(out) <- NULL
  return(out)
}
