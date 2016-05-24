
#############################################################################
plot.IRT.informationCurves <- function(x , curve_type="test" , ... ){

	theta <- x$theta[,1]
	args <- list(...)
	xlim <- base::range(theta)
	#************************************
    # collect arguments for plot
	if ( sum( names(args) == c("xlim") ) == 1 ){
			xlim <- args$xlim
						}	
	if ( sum( names(args) == c("xlim") ) == 0 ){
			xlim -> args$xlim
						}
	if ( sum( names(args) == "xlab") == 0 ){
			args$xlab <- "theta" 
						}
	if ( sum( names(args) == "ylab") == 0 ){
			args$ylab <- "y" 
						}		
					
					
	#************************************************************	
	#***** test information curve or se curve
	if ( curve_type %in% c("test" , "se" ) ){
		
		if ( curve_type == "se" ){	y <- x$se_curve }
		if ( curve_type == "test" ){	y <- x$test_info_curve }
        	    
		y1 <- y[ ( theta > xlim[1] ) & ( theta < xlim[2] ) ]
		if ( is.null(args$ylim) ){
			ylim <- c( 0 , max( y1 , na.rm=TRUE ) )
			args$ylim <- ylim
							}
		args$x <- theta
		args$y <- y
		args$type <- "l"	
        do.call( graphics::plot , args )		
	
						}
	#************************************************************
	#*********
	# add item information curves here ...

					}
#############################################################################					