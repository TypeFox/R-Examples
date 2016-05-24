
tetrachoric2 <- function( dat , method="Ol" ,  delta=.007 , maxit = 1000000 ,
	cor.smooth=TRUE , progress=TRUE){
	# process data
	dat <- as.matrix(dat)
	if (method == "Ol"){
		res <- polychoric2( dat=dat , cor.smooth=cor.smooth , maxiter=100)
						}		
	if ( method != "Ol" ){
	
		# calculate tau
		tau <- - stats::qnorm( colMeans(dat,na.rm=TRUE ) )
		dat.resp <- 1 - is.na(dat)
		dat[ dat.resp==0] <- 9
		I <- ncol(dat)
		# calculate frequencies
		dfr <- data.frame( "item1" = rep(1:I,I) , "item2" = rep(1:I, each=I ) )
		h1 <- t( ( dat==1 ) ) %*% ( dat==0 )
		dfr$f11 <- matrix( t( ( dat==1 ) ) %*% ( dat==1 ) , ncol=1 , byrow=TRUE ) + .5
		dfr$f10 <- matrix( t( ( dat==1 ) ) %*% ( dat==0 ) , ncol=1 , byrow=TRUE ) + .5
		dfr$f01 <- matrix( t( ( dat==0 ) ) %*% ( dat==1 ) , ncol=1 , byrow=TRUE ) + .5
		dfr$f00 <- matrix( t( ( dat==0 ) ) %*% ( dat==0 ) , ncol=1 , byrow=TRUE ) + .5

		# guessing
	#	if ( ! is.null(guess) ){
	#	   dfr0 <- dfr
	#	   dfr0$ftot <- dfr$f11 + dfr$f10 + dfr$f01 + dfr$f00
	#	   dfr$f00 <- dfr0$f00 / ( 1 - guess[dfr$item1 ] ) /  ( 1 - guess[dfr$item2 ] )
	#	   dfr$f01 <- ( (1 - guess[dfr$item2 ])* dfr0$f01 - guess[dfr$item2] * dfr0$f00 ) / 
	#						( 1 - guess[dfr$item1 ] ) /  ( 1 - guess[dfr$item2 ] )	
	#	   dfr$f10 <- ( (1 - guess[dfr$item1 ])* dfr0$f10 - guess[dfr$item1] * dfr0$f00 ) / 
	#						( 1 - guess[dfr$item1 ] ) /  ( 1 - guess[dfr$item2 ] )	
	#	   dfr$f11 <- dfr0$ftot - dfr$f00 - dfr$f01 - dfr$f10 - dfr$f11				
	#		}
			
		dfr$ftot <- dfr$f11 + dfr$f10 + dfr$f01 + dfr$f00
		dfr$p11 <- dfr$f11 / dfr$ftot
		dfr$pi1 <- ( dfr$f11 + dfr$f10 - 1 ) / dfr$ftot
		dfr$pi2 <- ( dfr$f11 + dfr$f01 - 1) / dfr$ftot
		# subdata of dfr
		dfr <- dfr[ dfr$item1 > dfr$item2 , ]
		dfr <- dfr[ dfr$ftot > 0 , ]	
		
		dfr$qi1 <- stats::qnorm( dfr$pi1)
		dfr$qi2 <- stats::qnorm( dfr$pi2)
		#******************
		# method of Bonett
		if ( method %in% c("Bo","Di") ){ 
			dfr$pmin <- ifelse( dfr$pi1 < dfr$pi2 , dfr$pi1 , dfr$pi2 )
			dfr$c <- ( 1 - abs( dfr$pi1 - dfr$pi2 ) / 5 - ( 0.5 - dfr$pmin)^2  ) / 2
			dfr$omega <- ( dfr$f00 * dfr$f11 ) / ( dfr$f01 * dfr$f10)
			dfr$r0 <- base::cos( pi / ( 1 + dfr$omega^( dfr$c )  ) )
						}
		#*****
		# method of Divgi
		if ( method == "Di"){	
			dfr2 <- as.matrix(dfr)
			numdiffparm <- .000001
			maxiter <- 100
	#		dfr$r0 <- tetrachoric2_rcpp_aux( dfr2 , numdiffparm , maxiter )	
			dfr$r0 <- .Call("tetrachoric2_rcpp_aux" ,
						   dfr2 , numdiffparm , maxiter , PACKAGE="sirt" )	
				}					
		#******************	
		# method of Tucker
		if ( method=="Tu"){
			
			# functions defined by Cengiz Zopluoglu   
		#	L <- function(r,h,k) {(1/(2*pi*sqrt(1-r^2)))*exp(-((h^2-2*h*k*r+k^2)/(2*(1-r^2))))}
			L <- function(r,h,k) {(1/(2*pi*sqrt(1-r^2)))*exp(-((h^2-2*h*k*r+k^2)/(2*(1-r^2))))}
			S <- function(x,r,h,k) { x-(L(r,h,k)*delta) }
			A0 <- dfr$A0 <- dfr$p11 - dfr$pi1 * dfr$pi2 
			dfr$r0 <- delta / 2
			dfr$iter <- 0

			dfr$iter <- 0
			dfr$conv <- 0
			ii <- 0

			vars <-  c("A0","r0","iter","conv")
			if (progress){ cat("    0 : " ) }
			while( ( mean( dfr$conv) < 1 ) & ( ii < maxit ) ){
				# iterations
				dfr0 <- dfr
				#***
				ii <- ii+1
				dfr$A0 <- dfr$A0 - delta * L( dfr$r0 , dfr$qi1 , dfr$qi2 )
				dfr$r0 <- dfr$r0 + delta
				dfr$iter <- ii
				ind <- which( dfr$A0 < 0 )
				if ( length(ind) > 0 ){ 
					h1 <- dfr$r0 - delta/2 + dfr$A0 / L( dfr0$r0 , dfr$qi1 , dfr$qi2 )
					dfr$r0[ind] <- h1[ind]
					dfr$conv[ind] <- 1
									}
				i2 <- which( dfr$conv==1 & dfr0$conv==1 )
				if (length(i2) > 0){    dfr[ i2 , vars] <- dfr0[ i2 , vars] }
				if (progress){
					if (ii %% 100 == 0 ){ cat(".") ; utils::flush.console() }
					if (ii %% 1000 == 0 ){ cat("\n" ,ii, ": ") }    
							}
					}
		#******************			
			cat("\n")
			}
		TC <- matrix(NA , I , I )
		diag(TC) <- 1
		TC[ as.matrix(dfr[ , c("item1","item2") ] ) ] <- dfr$r0
		TC[ as.matrix(dfr[ , c("item2","item1") ] ) ] <- dfr$r0
		if (cor.smooth){ 
			TC <- psych::cor.smooth(TC) 
					}
		rownames(TC) <- colnames(TC) <- colnames(dat)
		res <- list("tau"=tau , "rho" = TC )
		
		}   # method != "Ol"		
	return(res)
	}