"lines.cjs" <- 
function( x, what="n", animals=-1, occasions=-1, ... ){

if( what == "s" ){
	y <- x$s.hat
	if( animals[1] <= 0 ){
		animals <- 1:nrow(y)
	}
	occasion <- 1:nrow(y)
	if( occasions[1] <= 0 ){ 
		occasions <- 1:ncol(y)
	}
	occasions <- occasions[ occasions != ncol(y) ]
	y <- y[,occasions]
	occasion <- occasion[occasions]

	if( (length(animals) != 1) | (animals[1] > 0 ) ){
		y <- y[animals,]
	}

	matlines( occasion, t(y), ... )
} else {

	y <- x$n.hat
	occasion <- seq( along=y )
	if( occasions[1] <= 0 ){ 
		occasions <- 1:length(y)
	}
	occasions <- occasions[ occasions != 1 ]

	y <- y[occasions]
	occasion <- occasion[occasions]

	lines( occasion, y,  ...)

}

invisible(1)

}

