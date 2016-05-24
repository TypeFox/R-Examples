
##############################################
# residualize data
lsem.residualize <- function( data , moderator , moderator.grid ,
		lavmodel , h = 1.1 , residualize=TRUE , eps=1E-8 , verbose=TRUE ){	
    # lavaanify model
	lavaanstr <- lavaan::lavaanify( lavmodel  )
	vars <- unique( c( lavaanstr$rhs  , lavaanstr$lhs ) )
	vars <- intersect( colnames(data) , vars )
	data.mod <- data[ , moderator ]
	N <- length(data.mod)
	# select nearest neighbor in moderator group for calculating residuals
	G <- length(moderator.grid)
    modgrid_index <- rep(1,N)	
	for (gg in 2:G){
		modgrid_index <- ifelse( abs( data.mod - moderator.grid[ modgrid_index ] ) < 
					abs( data.mod - moderator.grid[ gg ] ) ,
					modgrid_index , gg )
		
			}
	# compute weights for every grid point gg
	weights <- matrix( NA , nrow=N , ncol=G )	
	sd.moderator <- stats::sd( data.mod , na.rm=TRUE)
	bw <- h * sd.moderator * N^(-1/5) 
	moderator.density <- stats::density( data.mod , from= min( moderator.grid) ,
				to = max(moderator.grid ) , n = G )$y		
	moderator.density <- data.frame( "moderator"=moderator.grid , 
				"wgt"=moderator.density / sum(moderator.density) )
	
	for (gg in 1:G){
		# gg <- 1
		xgg <- moderator.grid[gg]	
		wgt <- stats::dnorm( data.mod , mean = xgg , sd = bw ) / 
					stats::dnorm( xgg , mean=xgg , sd=bw )
		weights[,gg] <- ifelse( wgt < eps , eps , wgt )
				}
	
	dat2 <- data
	V <- length(vars)
	
	residualized_interceps <- matrix( 0 , nrow=G , ncol=V)
	colnames( residualized_interceps ) <- vars
	rownames( residualized_interceps ) <- round( moderator.grid , 3 )


	
	if ( residualize){
		if (verbose){ cat("** Residualize Data\n") ; utils::flush.console() }
		
		for (vv in 1:V){
		# vv <- 1
		var.vv <- vars[vv]	
		for (gg in 1:G){
			# gg <- 1
			x <- dat2[,moderator]
			mod <- stats::lm( data[ , var.vv ] ~  x + I( x^2 ) ,weights= weights[,gg]  )
			m1 <- stats::predict( mod , data.frame( x = moderator.grid[gg] ) )
			residualized_interceps[gg,vv] <- m1
			y <- stats::resid(mod)
			
			dat2[ , var.vv ] <- ifelse( modgrid_index == gg , y , dat2[ , var.vv ] )
						}	
					}
				}
				
	res <- list( "resid_vars" = vars , "data" = dat2 ,
			"weights_grid" = weights , "bw" = bw ,
			"moderator.density"=moderator.density ,			
			"sd.moderator"= sd.moderator , "G"= G , "N"=N ,
			"residualized_interceps" = residualized_interceps )
	return(res)
		}
###############################################		