findClosestLargerPowerOf2 <- function( x ){
	
	if( x > 0 ){
		
		out <- ceiling( log( x, 2 ) )
		return( out )
		
	} else{
		
		stop( "x must be a positive number." )
		invisible( NULL )
		
	}
	
}

###

normalizeToUnitInterval <- function( x ){
	
	## function to normalize the domain of an observed
	## curve to the unit interval
	
	range.x <- range( x )
	
	output <- ( x - range.x[1] ) / ( range.x[2] - range.x[1] )
	
	return( output )
	
}

###

projectCurveWavelets <- function( x, y, irreg.grid=FALSE, grid.length=NULL,
filter.number=10, family="DaubLeAsymm", bc="periodic", verbose=FALSE, ... ){
	
		## REQUIRES 'wavethresh' PACKAGE
	
		## normalize to unit interval
		x <- normalizeToUnitInterval( x )
		length.x <- length( x )
		
		## irregular grid?
		if( irreg.grid ){
			
			if( is.null( grid.length ) ){
				
				if( log( length.x, 2 ) %% 1 == 0 ){
					
					grid.length <- length.x
					
				} else{
					
					closest.power <- findClosestLargerPowerOf2( length.x )
					grid.length <- 2^closest.power
								
				}
					
			} else if( grid.length <= 0 || ( log( grid.length, 2 ) %% 1 != 0 ) ){

				stop( "grid.length must be a positive power of 2.")				
				
			}
				
			## make regular grid
			grid <- makegrid( x, y, gridn=grid.length )
			x <- grid$gridt
			y <- grid$gridy
				
			## wavelet transform
			y.wd <- irregwd( gd=grid, filter.number=filter.number, family=family, 
			bc=bc, verbose=verbose )					
				
		} else if( is.null( grid.length ) ){
			
			if( log( length.x, 2 ) %% 1 == 0  ){
				
				## wavelet transform
				y.wd <- wd( y, filter.number=filter.number, family=family, 
				bc=bc, verbose=verbose )
				
			} else{
				
				closest.power <- findClosestLargerPowerOf2( length.x )
				grid.length <- 2^closest.power

				## make regular grid
				grid <- makegrid( x, y, gridn=grid.length )
				x <- grid$gridt
				y <- grid$gridy

				## wavelet transform
				y.wd <- irregwd( gd=grid, filter.number=filter.number, family=family, 
				bc=bc, verbose=verbose )					
				
			}
		
		} else if( grid.length <= 0 || ( log( grid.length, 2 ) %% 1 != 0 ) ){
			
				stop( "grid.length must be a positive power of 2.")							
			
		} else{
			
			## make regular grid
			grid <- makegrid( x, y, gridn=grid.length )
			x <- grid$gridt
			y <- grid$gridy			
			
			## wavelet transform
			y.wd <- irregwd( gd=grid, filter.number=filter.number, family=family, 
			bc=bc, verbose=verbose )					
			
		}
						
		## thresholding
		y.wdT <- threshold( y.wd, ... )
		
		## extract coefficients
		y.coeffs <- y.wdT$D
		
		## smooth curve
		y.wave <- wr( y.wdT )
		
		output <- list( coefficients=y.coeffs, y.wdT=y.wdT, y.wavelet=y.wave,
		x.grid=x )
		
		return( output )
		
}