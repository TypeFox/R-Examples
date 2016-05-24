


#########################################################
# create A matrix
tamaanify.create.A <- function( res ){
    resp <- res$resp
	lavpartable <- res$lavpartable
	#***********************************
	# create A matrix	
	A <- designMatrices( resp=resp )
	A <- A$A
	items0 <- dimnames(A)[[1]] <- res$items$item	
	
	K2 <- dim(A)[2]
	if (K2==2){
	   dimnames(A)[[3]] <- paste0( dimnames(A)[[3]] , "_Cat1")
				}

				
	#*********************************************
	# loop over items for smoothed nominal response models
	smooth.nrm <- FALSE
	itemtable <- res$items	
	items.ind <- grep( "," , paste(itemtable$itemtype) , fixed=TRUE )
	for (ii in items.ind ){	
		# ii <- ind[1]
		smooth.nrm <- TRUE
		itemtype.ii <- paste( itemtable[ ii , "itemtype"] )
		itemtype.ii <- strsplit( itemtype.ii , split="," , fixed=TRUE )[[1]]
		maxK.ii <- itemtable[ii,"ncat"] - 1
		if ( itemtype.ii[1] %in% c("GPCM" , "PCM" )){
			n_order <- as.numeric( itemtype.ii[2] )
						}
		if ( ! ( itemtype.ii[1] %in% c("GPCM" , "PCM" ) ) ){
			stop(paste0("Only item types 'GPCM' or 'PCM' are allowed \n  for ",
					"smoothing item intercepts.") )
						}
		# extend it to nominal response model!!
		item.ii <- paste(itemtable[ii , "item" ])
		ind1 <- which( dimnames(A)[[3]] == paste0( item.ii , "_Cat1") )
		A[ ii , 1 + 0:maxK.ii , ind1 ] <- - ( 0:maxK.ii )
		dimnames(A)[[3]][ind1] <- paste0( item.ii , "_lin" )
		
		# further fourier terms
		if ( n_order > 1){
			for (ff in 2:n_order){
		#		ff <- 2
				ind1 <- which( dimnames(A)[[3]] == paste0( item.ii , "_Cat",ff) )
				A[ ii , 1 + 0:maxK.ii , ind1 ] <- 
						- sin( 3.141593 *( 0:maxK.ii ) * (ff - 1 ) / maxK.ii )
				dimnames(A)[[3]][ind1] <- paste0( item.ii , "_four" , ff-1)
							}
						}
		if ( n_order < maxK.ii ){
			for (ff in (n_order+1):maxK.ii){
				ind1 <- which( dimnames(A)[[3]] == paste0( item.ii , "_Cat",ff) )
				var.ii <- dimnames(A)[[3]][ind1]
				A <- A[ , , -c(ind1) ]
							}	
						}
				}
			

	#**********************************		
	#*********************************************
	#****** xsi parameter fixings	
	xsi.fixed <- NULL
	maxK <- max( res$items$ncat )-1
	lavpartable <- lavpartable[ lavpartable$user != -99 , ]
	# handle xsi parameter fixings
	for (hh in 1:maxK){
		# hh <- 1
		vv <- paste0( "t" , hh )
		ind <- which( paste0( lavpartable$op , lavpartable$rhs ) == 
					paste0( "|" , vv )  )
		ind1 <- which( lavpartable$free == 0 ) 
		ind <- intersect( ind , ind1 )
		if ( length(ind) > 0 ){		
			lv1 <- lavpartable[ ind , ]
			N1 <- nrow(lv1)
			for (zz in 1:N1){
				# zz <- 1
				lv1.zz <- lv1[zz,]
				i1 <- which( items0 %in% paste( lv1.zz$lhs ) )
				Azz <- A[ i1 , hh+1 , ]
				i2 <- which( Azz != 0 )
				xsi.zz <- cbind( i2 , - lv1.zz$ustart )
				xsi.fixed <- rbind( xsi.fixed , xsi.zz )
						}
					}
				}
			
				
	#******************************
	#******************************
	# xsi equality constraints
	lavpartable <- res$lavpartable
	lavpartable0 <- lavpartable
	thresh <- paste0( "t" , 1:maxK )
	ind1 <- which( lavpartable0$op == "|" & ( lavpartable0$rhs) %in% thresh )
	ind2 <- which( paste(lavpartable0$label) != "" )
	ind <- intersect(ind1,ind2)
	lavpartable0 <- lavpartable[ind,]
	lavpartable0$label2a <- paste0( lavpartable0$lhs , "_Cat" ,
					 substring( lavpartable0$rhs , 2 ) )
	
	labs <- unique( paste0( lavpartable0$label ))	
	NL <- length(labs)

	for (ll in 1:NL){
		# ll <- 1	
		labs.ll <- labs[ll]
		lav.ll <- lavpartable0[ paste0(lavpartable0$label) == labs[ll] , ]
		t11 <- as.numeric( substring(lav.ll[ 1 , "rhs" ],2)) +1 
		i11 <- paste0(lav.ll$lhs)
		if ( length(i11) > 1 ){
			A00 <- A
			A10 <- A[ items0 %in% i11[1] , t11 , ] 
			i10 <- which( A10 != 0 )
			A00[ items0 %in% i11[-1] , t11 , i10 ] <- A10[ i10 ]
			A11 <- A[ items0 %in% i11[-1] , t11 , , drop=FALSE ] 
			A11 <- colSums(abs(A11))
			A <- A00	
			A <- A[ , , - which(A11>0) ]		
				}


		}
		ind <- match( paste(dimnames(A)[[3]]) , paste( lavpartable0$label2a) )
		if ( ! smooth.nrm ){
			dimnames(A)[[3]] <- paste(lavpartable0[ ind , "label" ])
						}
						
	#***********************
	# model constraint thresholds	
	
	mdfr <- res$MODELCONSTRAINT.dfr

	if ( ! is.null(mdfr) ){
		mdfr <- mdfr[ grep( "|t" , mdfr$fullsyn , fixed=TRUE) , ]
		lav1 <- res$lavpartable
		lav1 <- lav1[ grep( "|t" , lav1$fullsyn , fixed=TRUE) , ]
		if ( nrow(mdfr) > 0 ){
			items0 <- colnames(res$resp)
			lav1$itemno <- match( lav1$lhs , items0 )		
			lav1 <- lav1[ order(lav1$itemno) , ]		
			syn2 <- mdfr$syn
			# add parameters		
			dfr <- tamaanify.grep.linequations( syn2 )

			lav1a <- lav1[ match( dfr$lhsparm , paste(lav1$label ) ) , ]
			dfr <- cbind( lav1a[ , c("fullsyn" , "lhs" , "op" , "rhs" ) ] , dfr )
			dfr$user <- 1			
			dfr$ustart <- NA
			ind <- which( ! paste( lav1$label ) %in% dfr$lhsparm )
			if (length(ind) > 0 ){
				lav1 <- lav1[ ind , ]
				lav1$terms <- lav1$rhsparm <- lav1$lhsparm <- lav1$label
				lav1$fac <- 1
				lav1$parm <- lav1$label		
				lav1 <- lav1[ , colnames(dfr) ]
				dfr <- rbind( dfr , lav1 )					
						}
			rownames(dfr) <- NULL	
			res$modelconstraint.thresh <- dfr	
			#*** create A matrix according model constraints
			dfr <- dfr[ dfr$user > 0 , ]
			
			parms <- unique( paste(dfr$parm))
			Nparm <- length(parms)
			I <- ncol(res$resp)
			maxcat <- res$maxcat	
			A <- array( 0 , dim=c(I, maxcat+1 , Nparm) )
			dimnames(A)[[1]] <- colnames(res$resp)
			dimnames(A)[[2]] <- paste0("Cat" , 1:(maxcat+1))
			dimnames(A)[[3]] <- parms		
			ND <- nrow(dfr)
			for (dd in 1:ND){
				# dd <- 1
				dfr.dd <- dfr[dd,]
				idd <- as.numeric( substring( dfr.dd$rhs , 2 ) ) + 1
				A[ dfr.dd$lhs , idd , dfr.dd$parm ] <- - dfr.dd$fac
							}
							}
					}
	#**************************************
	#*** output
	res$A <- A
	res$xsi.fixed <- xsi.fixed	
	return(res)
	}
#########################################################
