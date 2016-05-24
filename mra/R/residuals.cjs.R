residuals.cjs <- function( object, type="pearson", ... ){

if( type == "deviance" ){
    hists <- object$histories
    hists[ hists >=2 ] <- 1  # Don't need '=' of '>=' because only 0,1, and 2 present, but included it just to be sure we only have 0's and 1's

    if( !is.null(object$fitted) ){
    	cellmeans <- object$fitted
    } else {
    	cellmeans <- predict.cjs( object )
    }

	resids <- -sqrt(2*abs( log(1-cellmeans) ))   # these are for hists == 0
	resids[ hists == 1 ] <- sqrt(2*abs(log(cellmeans[ hists == 1 ])))  # these are for hists == 1
} else {
    if( !is.null(object$residuals) ){
	   resids <- object$residuals 
    } else {
        hists <- object$histories
        hists[ hists >=2 ] <- 1  
    
        if( !is.null(object$fitted) ){
        	cellmeans <- object$fitted
        } else {
        	cellmeans <- predict.cjs( object )
        }
    
    	resids <- (hists - cellmeans) / sqrt(cellmeans*(1-cellmeans))
    }
}

resids
}

