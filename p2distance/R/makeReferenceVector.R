makeReferenceVector <-
function (X, reference_vector_function = min){
	
		vRef <- as.matrix(apply(X, MARGIN=2, reference_vector_function))	# Calcula el minimo de cada variable contenida en la matrix X
		vRef <- t(vRef)	# Transpone lo anterior 		
		
		return(vRef)			
}

