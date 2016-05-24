summary.irates <-
function( object, 
							... 
							){

  if(class(object) != "irates"){stop("Object needs to be of class irates")}
  cat("\n Baseline covariate distribution \n")
  basel = as.data.frame(object$N, row.names = object$covar.lab)
  names(basel) = "individuals"
  print(basel, quote=FALSE, right=FALSE, ...)
  
  cat("\n\n")

  print(x = object, 
		event.code = object$event.code,
		covar.code = object$covar.code, 
		full.sample = TRUE,
		display.digits = 4 
		)
							
}

