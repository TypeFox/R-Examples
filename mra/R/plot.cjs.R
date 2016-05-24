"plot.cjs" <- 
function( x, type="n", ci=TRUE, smooth=TRUE, occasions=-1, animals=-1,
	smubass=5, ... ){

	if( type == "s" ){
		#  Plot survival
		y <- x$s.hat
		occasion <- 1:ncol(y)
		if( any(occasions <= 0) ){ 
			occasions <- 1:ncol(y)
		}
		occasions <- occasions[ occasions != ncol(y) ]
		if( is.null(dimnames(y)[[2]]) ){
			nms <- 1:ncol(y)
		} else {
			nms <- dimnames(y)[[2]]
		}
		y <- y[,occasions]
		occasion <- occasion[occasions]
		nms   <- nms[occasions]
		if( (length(animals) != 1) | (animals[1] > 0 ) ){
			y <- y[animals,]
		}


		survival <- y
		if( is.matrix(y) ){
			if( nrow(y) > 1 ){ 
				survival <- t( y )
			}
		}


		matplot(occasion, survival, type="l", xaxt="n", xlab="Period", 
		ylab="Survival", ...)
		axis(1, at=occasion, labels=nms )
	
		ans <- survival
	} else {
        ans <- plot.nhat( x, ci, smooth, occasions, smubass, ... )
	}


	invisible( ans )
}

