

######################################################################
MHprop_hrm <- function( MHprop , b , a , phi , theta , iter , burnin ){

	psi <- phi
	MHprop <- list()
	
	# logical indicating whether standard deviation of the
	# proposal distribution is calculated according to
	# Browne and Draper (2006)
	MHprop$refresh_formula <- TRUE
	
	MHprop$SD <- list( "b" = .4+0*b , "a" = .4 +0*a , "phi" = .2 + 0*phi , "psi" = .4 + 0*psi ,
			"theta" = .5 + 0*theta ) 
	MHprop$accept <- list("b" = 0 + 0*b , "a" = 0 + 0*a , "phi" = 0+0*phi ,
						"psi" = 0 + 0*psi , "theta" = 0 + 0*theta  )

	MHprop$refresh_count$b <- 0
	MHprop$refresh_count$a <- 0
	MHprop$refresh_count$phi <- 0
	MHprop$refresh_count$psi <- 0
	MHprop$refresh_count$theta <- 0

	RI <- 25
	MHprop$refresh_iter$b <- RI
	MHprop$refresh_iter$a <- RI
	MHprop$refresh_iter$phi <- RI
	MHprop$refresh_iter$psi <- RI
	MHprop$refresh_iter$theta <- RI

	MHprop$refresh_SDchange$b <- .02
	MHprop$refresh_SDchange$a <- .02
	MHprop$refresh_SDchange$phi <- .02
	MHprop$refresh_SDchange$psi <- .02
	MHprop$refresh_SDchange$theta <- .05
	
	
	vars <- c("b" ,"a" , "phi" , "psi" , "theta") 
	
	# compute iterations for which MH updatings must be computed
	refresh_iters <- sort( unique( unlist( MHprop$refresh_iter ) ) )
	RI <- length(refresh_iters)
	v1 <- NULL
	for (rr in 1:RI){
		 l1 <- ( 1:iter %% refresh_iters[rr] ) == 0 
		 l1 <- (1:iter)[l1]
	     v1 <- c( v1 , l1 ) 
					}
	v2 <- sort( unique(v1) )
	v2 <- v2[ v2 <= burnin ]
	MHprop$ITER_refreshing <- v2
	
	
	# refreshing variables
	MHprop$VARS_refreshing <- c("b" , "phi" , "psi" , "theta") 
	
	# boundaries for acceptance rates
	MHprop$accept_bounds <- c( .4 , .6 )
	
	return(MHprop)
	}
######################################################################	