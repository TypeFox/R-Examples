computeLVConstraints <- function(results,x_axis=1,y_axis=2){
		if( ("lx" %in% names(results)) && ("ly" %in% names(results)) ){
			return(minmaxHelper(results$lx,results$ly,axis1=x_axis,axis2=y_axis,findBounds=FALSE))
		}else{
			stop("Results are improperly formatted. computeLVConstraints must stop.")
		}	
}