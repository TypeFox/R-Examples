
				
################################################################
# collect item parameters
.gdm.collect.itempars <- function( data , K , D , b , a , TD , thetaDes , irtmodel ,
			se.b , se.a , data0){	
		# collect item parameters	
		item <- data.frame( "item" = colnames(data0) ,
					"N" = colSums(1-is.na(data0) ) )
		item$M <- colMeans(data0 , na.rm=T)
		# b parameters
		se.b[ b < -9999  ] <- NA		
		b[ b < -9999  ] <- NA
		se.a[ a < -9999  ] <- NA		
		a[ a < -9999  ] <- NA		
		
		for (kk in 1:K){
		item[ , paste0( "b.Cat" , kk) ] <- b[,kk]
						}						
		for (dd in 1:TD){
			if ( irtmodel %in% c("1PL" , "2PL") ){		
				item[ , paste0( "a." , colnames(thetaDes)[dd] ) ] <- a[,dd,1]
								}
			if ( irtmodel %in% c("2PLcat") ){		
			  for (kk in 1:K){ 
					item[ , paste0( "a." , colnames(thetaDes)[dd] ,
								".Cat", kk ) ] <- a[,dd,kk]								
									}
		se.a[ a == -99999  ] <- NA		
		a[ a == -99999  ] <- NA
								}
						}
		res <- list("item" = item , "b" = b , "se.b" = se.b , "a" = a )
		return(res)
			}
#######################################################################
# moments of distribution
.gdm.calc.distributionmoments <- function( D , G , pi.k , theta.k ){
	mean.trait <- sd.trait <- skewness.trait <- matrix( 0 , nrow=D , ncol=G )			
	for (dd in 1:D){
		for (gg in 1:G){
		mean.trait[dd,gg] <- sum( theta.k[,dd] * pi.k[ , gg ] )
		sd.trait[dd,gg] <- sqrt( sum( theta.k[,dd]^2 * pi.k[ , gg ] ) - mean.trait[dd,gg]^2 )
		skewness.trait[dd,gg] <- sum( ( theta.k[,dd] - mean.trait[dd,gg] )^3 * pi.k[ , gg ] ) /
						sd.trait[dd,gg]^3
					}
			}
	rownames(skewness.trait) <- rownames(sd.trait) <- 
				rownames(mean.trait) <- colnames(theta.k)
	colnames(skewness.trait) <- colnames(sd.trait) <- 
				colnames(mean.trait) <- paste0("Group",1:G)
	#*****
	# correlation matrices
	correlation.trait <- as.list(1:G)
	names(correlation.trait) <- colnames(mean.trait)
	for (gg in 1:G){
		# gg <- 1
		mean.gg <- rep(0,D)
		Sigma.gg <- diag(0,D)	
		for (dd in 1:D){
			mean.gg[dd] <- sum( pi.k[,gg] * theta.k[,dd] )
				}
		for (dd1 in 1:D){
			for (dd2 in dd1:D){
#			for (dd2 in 1:D){
		#		dd1 <- 1 ; 	dd2 <- 1
				Sigma.gg[dd1,dd2] <- sum( pi.k[,gg] * (theta.k[,dd1] - mean.gg[dd1] )*(theta.k[,dd2] - mean.gg[dd2] ) ) 
#				Sigma.gg[dd1,dd2] <- Sigma.gg[dd1,dd2] - mean.gg[dd1] * mean.gg[dd2]
				Sigma.gg[dd2,dd1] <- Sigma.gg[dd1,dd2]
									}
						}
		rownames(Sigma.gg) <- colnames(Sigma.gg) <- rownames(mean.trait)	
		correlation.trait[[gg]] <- stats::cov2cor(Sigma.gg + diag(10^(-5),D) )
		
					}	
	res <- list( "mean.trait"=mean.trait , "sd.trait" = sd.trait , 
				"skewness.trait" = skewness.trait , "correlation.trait"=correlation.trait)
    return(res)				
		}

#############################################################
# calculation of information criteria and number of parameters
.gdm.calc.ic <- function( dev , dat , G ,  skillspace , irtmodel , 
			K,D,TD,I,b.constraint,a.constraint , mean.constraint ,
			Sigma.constraint , delta.designmatrix , standardized.latent ,
			data0 , centerslopes , TP , centerintercepts ){
    ic <- list( "deviance" = dev , "n" = nrow(data0) )
	ic$traitpars <- 0
	ic$itempars <- 0	
	#******
	# Up to now this works in one dimension
	# trait parameters: normal skillspace
	if ( skillspace == "normal" ){
		if (irtmodel=="1PL" & ( D==1 )){
			ic$traitpars <- 2*(G-1)	+ 1
							}
														
		if ( ( irtmodel %in% c("2PL","2PLcat") ) & (D==1) ){
			ic$traitpars <- 2*(G-1)
							}	
		if (D > 1 ){
			ic$traitpars <- 2 * D*G + D*(D-1)/2*G
			if ( ! is.null(mean.constraint) ){ 
					ic$traitpars <- ic$traitpars - nrow(mean.constraint)			
									}
			if ( ! is.null(Sigma.constraint) ){ 
					ic$traitpars <- ic$traitpars - nrow(Sigma.constraint)			
									}
#			if (standardized.latent ){ ic$traitpars <- ic$traitpars - D }									
							}														
						}	# end normal
	#******
	# trait parameters: loglinear skillspace
	if ( skillspace == "loglinear" ){
			ic$traitpars <- G*(ncol(delta.designmatrix) - 1)
						}
	if ( skillspace == "full" ){
			ic$traitpars <- G*(TP-1)
						}	
	if ( skillspace == "est" ){
			ic$traitpars <- G*(TP-1) + TP*TD
						}							
	#************************************************
	# item parameters b
	ic$itempars.b <- I*K
	if ( ! is.null(b.constraint)){
			ic$itempars.b <- ic$itempars.b - nrow(b.constraint)
							}	
	#************************************************
	# item parameters a
	ic$itempars.a <- 0
	if ( irtmodel == "2PL"){ 
		ic$itempars.a <- I*TD
		if ( ! is.null(a.constraint)){
				a.constraint2 <- a.constraint[ a.constraint[,3] == 1 , ]
				ic$itempars.a <- ic$itempars.a - nrow(a.constraint2)
								   }	
						}
	ic$centeredintercepts <- (centerintercepts)*D						
	ic$centeredslopes <- (centerslopes)*D
	if ( irtmodel == "2PLcat"){ 
		ic$itempars.a <- I*TD*K
		if ( ! is.null(a.constraint)){
				ic$itempars.a <- ic$itempars.a - nrow(a.constraint)
								   }	
						}
	#***********************************************
	# information criteria
	ic$itempars <- ic$itempars.a + ic$itempars.b - ic$centeredintercepts - ic$centeredslopes
	ic$np <- ic$itempars + ic$traitpars	
#	ic$n <- n # number of persons
	# AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC (conistent AIC)
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )				
    return(ic)
		}
###################################################################
		
###########################################
# person parameter estimates
.gdm.person.parameters <- function( data , D , theta.k ,
	p.xi.aj , p.aj.xi , weights ){
	#**************************
	person <- data.frame("case" = 1:(nrow(data)) , "M" = rowMeans( data , na.rm=T)  )
    EAP.rel <- rep(0,D)
    names(EAP.rel) <- colnames(theta.k)
    nstudl <- rep(1,nrow(data))
	doeap <- TRUE
	if ( is.list( p.aj.xi)){	
		p.aj.xi <- p.aj.xi[[1]] 
		nstudl <- rep(1,nrow(p.aj.xi ) )
		weights <- weights[,1]  
		weights <- weights[ weights > 0 ]
		doeap <- FALSE
			}
	if (doeap ){
	for (dd in 1:D){ #dd <- 1
		dd1 <- colnames(theta.k)[dd]
		person$EAP <- rowSums( p.aj.xi * outer( nstudl , theta.k[,dd] ) )
		person$SE.EAP <- sqrt(rowSums( p.aj.xi * outer( nstudl , theta.k[,dd]^2 ) ) - person$EAP^2)	
		EAP.variance <- stats::weighted.mean( person$EAP^2 , weights ) - ( stats::weighted.mean( person$EAP , weights ) )^2
		EAP.error <- stats::weighted.mean( person$SE.EAP^2 , weights )
		EAP.rel[dd] <- EAP.variance / ( EAP.variance + EAP.error )	
		colnames(person)[ which( colnames(person) == "EAP" ) ] <- paste("EAP." , dd1 , sep="")
		colnames(person)[ which( colnames(person) == "SE.EAP" ) ] <- paste("SE.EAP." , dd1 , sep="")				
		}

	# MLE	
	mle.est <- theta.k[ max.col( p.xi.aj ) , , drop=FALSE]
	colnames(mle.est) <- paste0( "MLE." , names(EAP.rel))
	person <- cbind( person, mle.est )
	# MAP
	mle.est <- theta.k[ max.col( p.aj.xi ) , , drop=FALSE]
	colnames(mle.est) <- paste0( "MAP." , names(EAP.rel))
	person <- cbind( person, mle.est )
			}
	# results	
	res <- list( "person" = person , "EAP.rel" = EAP.rel )
	return(res)
		}
###########################################################################	