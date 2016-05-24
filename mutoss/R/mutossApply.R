mutoss.apply <- function(mutossObj, f, label = deparse(substitute(f)), recordHistory = TRUE, ...) {
  
  
	params <- list()	
	for (param in names(formals(f))) {					#runs over all parameters of f
		if (param %in% slotNames(mutossObj)) {			#checks if parameter name corresponds to a mutoss slot	
			paramTail <- list(slot(mutossObj, param))	#extracts values from appropriate slot
			names(paramTail) <- param 					#attaches names for parameter values
			params <- c(params, paramTail)				#concatenates parameters with and without data from objects
		}
	}	
	result <- eval(as.call(c(f,params,...)))			#evaluates function call with extracted parameters	
	for (param in names(result)) {						
		if (param %in% slotNames(mutossObj)) {
                        value <- result[param][[1]]								
			attr(value, "method.name") <- label			#attaches attributes to 
			slot(mutossObj, param) <- value				#writes result in corresponding mutoss slots
		}
	}

        if ( recordHistory ) {
        
        mutossObj@commandHistory = c( mutossObj@commandHistory,
          paste(format( match.call( call = sys.call( sys.parent(1) ))), collapse='')) #write command into history
      }

	return(mutossObj)
}

