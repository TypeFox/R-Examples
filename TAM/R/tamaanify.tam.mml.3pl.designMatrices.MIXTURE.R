

######################################
# MIXTURE
tamaanify.tam.mml.3pl.designMatrices.MIXTURE <- function( res ){
	anlist <- res$ANALYSIS.list
	items <- colnames(res$resp)
	I <- length(items)
	itemtable <- res$items

    m1 <- mean( itemtable$itemtype %in% c("PCM","Rasch") )
	raschtype <- if ( m1 == 1 ){ TRUE } else { FALSE }
# Revalpr("raschtype")

	gammaslope.fixed <- NULL

		Q <- res$Q
		A <- res$A		
		D <- ncol(Q)   # number of dimensions

	  ncl <- anlist$NCLASSES
	  if (D==1){
		theta.grid <- seq( -6 , 6 , len=15 )
				}
	  if (D %in% c(2,3) ){
		theta.grid <- seq( -6 , 6 , len=10 )
				}		
	  if (D >= 4 ){
		theta.grid <- seq( -6 , 6 , len=5 )
				}		
				
	  # theta.grid <- seq( -4 , 4 , len=13 )
	  
	  ndim <- ncol(Q)
	  nodes <- theta.grid
	  theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , 
			ncol = ndim ) ) ) )
	  colnames(theta) <- colnames(Q)
	  TG <- nrow(theta)
	  TP <- TG * ncl
	  
	  res$theta.k <- diag(TP)
	  res$theta_MIXTURE <- theta	
	  
	   # design matrix E 
	   # E ( item , category , skill class , parameter )	

		
		#***** number of parameters to be estimated
		lavpartable <- res$lavpartable

		# estimated loadings in Q
		qloads <- sum( Q != 0 )
		Nparm <- dim(A)[3] * ncl + qloads * ncl

#		Nparm <- TP*D + sum( itemtable$ncat - 1 ) + 1		
		maxK <- res$maxcat
		
		# data frame containing estimated parameters
		parmsA <- as.vector(dimnames(A)[[3]])
		
		lav1 <- lavpartable[ lavpartable$op == "=~" , ]
		
		#--
		# intercepts
		p1 <- rep(parmsA , ncl)
		LA <- length(parmsA)
		dfr1 <- data.frame( "parm" =  paste0( p1 , "_Cl" , rep( 1:ncl , each=LA) ) ,
					"Class" = rep(1:ncl, each=LA) )
		dfr1$Cat <- lapply( strsplit( paste(dfr1$parm) , split="_" ) , FUN = function(pp){
						pp[2] } )
		dfr1$Cat <- gsub( "Cat" , "" , dfr1$Cat )				
		dfr1$dim <- NA
		dfr1$int <- 1
		dfr1$slo <- 0
		dfr1$fixed <- 0		
		dfr1$val <- NA
		dfr1$Aindex <- rep( seq(1,LA) , ncl )
		dfr <- dfr1
				
		#--
		# slopes
		p1 <- paste( lav1$label )
		dfr1 <- data.frame( "parm" =  paste0( p1 , "_Cl" , rep( 1:ncl , each=LA) ) ,
					"Class" = rep(1:ncl, each=LA) )
		dfr1$Class <- rep(1:ncl, each=LA)
		dfr1$Cat <- NA				
		dfr1$dim <- rep( paste(lav1$lhs) , ncl )
		dfr1$int <- 0
		dfr1$slo <- 1
		dfr1$fixed <- 1*( rep( lav1$free == 0 , ncl ) )
		dfr1$val <- rep( lav1$ustart , ncl )
		dfr1$Aindex <- NA
		dfr <- rbind( dfr , dfr1 )
		dfr$item <- lapply( strsplit( paste(dfr1$parm) , split="_" ) , FUN = function(pp){
						pp[1] } )
		dfr$index <- 1:(nrow(dfr))
		itempartable <- res$items
		item1 <- itempartable[ itempartable$itemtype %in% c("Rasch","PCM") , ]		
		ind1 <- which( paste(dfr$item) %in% paste(item1$item)  )
		ind2 <- which( dfr$slo == 1 )
		ind <- intersect( ind1 , ind2 )
		if ( length(ind) > 0 ){	
				dfr$fixed[ ind ] <- 1
							}
		res$itempartable_MIXTURE <- dfr
		
		#*******************
		# define E parameter design matrix
		
	   # design matrix E 
	   # E ( item , category , skill class , parameter )		
	   Nparm <- nrow(dfr)	
	   E <- array( 0 , dim= c(I , maxK+1 , TP , Nparm ) )	
	   dimnames(E)[[1]] <- colnames(res$resp)
	   p1 <- rep( paste0("theta" , 1:TG ) , ncl )
	   p1 <- paste0( p1 , "_Cl" , rep( 1:ncl , each = TG) )
	   dimnames(E)[[3]] <- p1 
	   dimnames(E)[[4]] <- paste(dfr$parm)
	   
	   #***
	   # create design matrix
	   
	   #- intercepts
	   dfr11 <- dfr[ dfr$int == 1 , ]
	   for (cl in 1:ncl){
		   # cl <- 2
		   # cl_temp <- paste0( "_Cl" , cl )
		   dfr11c <- dfr11[ dfr11$Class == cl , ]
		   for (hh in 1:(maxK+1) ){
			   # hh <- 1
			   for (tt in c( 1:TG + TG*(cl-1) ) ){
			       E[ , hh , tt, dfr11c$index ] <- 
						A[ ,hh , dfr11c$Aindex ]
										}
									}
						}
						
		#-- slopes

	   dfr11 <- dfr[ dfr$slo == 1 , ]
	   for (cl in 1:ncl){
	#		    cl <- 2
			   # cl_temp <- paste0( "_Cl" , cl )
			dfr11c <- dfr11[ dfr11$Class == cl , ]		
			NC <- nrow(dfr11c)
			
			for (nn in 1:NC){
			for (hh in 1:maxK){		
			dd <- paste( dfr11c$dim[nn] )
			for (tt in 1:TG ){
				tt1 <- tt + (cl-1)*TG
				E[ paste(dfr11c$item[nn]) , hh + 1 , tt1, dfr11c$index[nn] ] <- 
					E[ paste(dfr11c$item[nn]) , hh + 1 , tt1, dfr11c$index[nn] ] +
						dfr11c$val[nn] * hh * theta[ tt , dd] 
											}
								}
							}
					}
			#-- fixed gammaslope parameters
			dfr11 <- dfr[ dfr$slo == 1 , ]
			dfr11 <- dfr11[ dfr11$fixed == 1 , ]
			p1 <- dfr11[ , c( "index" , "val" ) ]
			colnames(p1) <- NULL
			gammaslope.fixed <- rbind( gammaslope.fixed , p1 )
			
			#*** constraints intercepts
			
			# zero sum constraint on item difficulties
			dimE4 <- dim(E)[4]
			
			ncolV <- ncl
			if ( ! raschtype ){
			   nQ <- ncol(Q)
			   ncolV <- ncolV + ncl*nQ
							}
			
			gammaslope.constr.V <- matrix( 0 , nrow=dimE4 , ncol=ncolV )
			rownames(gammaslope.constr.V) <- dimnames(E)[[4]]
			gammaslope.constr.c <- rep(0, ncl)
			if ( ! raschtype ){
				gammaslope.constr.c <- c(gammaslope.constr.c  , rep(1, ncl) )
							}
			

			
			dfr11 <- dfr[ dfr$int == 1 , ]			
			for (cl in 1:ncl){
				# cl <- 1
				dfr11c <- dfr11[ dfr11$Class == cl , ]
				gammaslope.constr.V[ dfr11c$index , cl ] <- 1
							}
            if ( ! raschtype ){
			    Qnames <- colnames(Q)
				dfr11 <- dfr[ dfr$slo == 1 , ]			
				ii <- ncl+1
				for (cl in 1:ncl){
					# cl <- 1
					dfr11c <- dfr11[ dfr11$Class == cl , ]
					for ( qname in Qnames ){
					   dfr11c1 <- dfr11c[ dfr11c$dim == qname , ]
   					   gammaslope.constr.V[ dfr11c1$index , ii ] <- 1 / nrow(dfr11c1)
					   ii <- ii+1
					                         }
								}
							}
							
	   res$raschtype <- raschtype		
	   res$gammaslope.fixed <- gammaslope.fixed
	   res$gammaslope.constr.V <- gammaslope.constr.V
	   res$gammaslope.constr.c <- gammaslope.constr.c
	   res$E <- E	
	   
	 return(res)
		}
