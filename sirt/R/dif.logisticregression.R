 


##NS export(dif.logistic.regression)
#---------------------------------------------------------------------------------------##
# This function performs itemwise DIF analysis by using logistic regression methods     ##
# uniform and nonuniform DIF                                                            ##
dif.logistic.regression <- function( dat , group , score ,
		quant=1.645){
    # INPUT:
    # dat       ... data frame (must only include item responses)
    # group     ... group identifier (this has to be a dummy variable)
    # score     ... matching criterion
    
    I <- ncol(dat)
    matr <- NULL
	cat("Items ")
    for (ii in 1:I){
     # ii <- 6
        dat.ii <- na.omit(data.frame( "y" = dat[,ii] , "score" = score , "group" = group ))
        mod1 <- stats::glm( y  ~ score , data = dat.ii , family="binomial")
        mod2 <- stats::glm( y  ~ score + group , data = dat.ii , family="binomial")
        mod3 <- stats::glm( y  ~ score + group + score*group , data = dat.ii , family="binomial")

		h1 <- data.frame( "item" = colnames(dat)[ii] , 
				"N" = sum( 1- is.na( dat[,ii] ) , na.rm=T) )
		h1$R <- min(group)
		h1$F <- max(group)
		h1$nR <- sum(  ( 1- is.na( dat[,ii] ) )* (1-group) , na.rm=T)
		h1$nF <- sum(  ( 1- is.na( dat[,ii] ) )* (group) , na.rm=T)		
		h1$p <- mean(  dat[,ii], na.rm=T) 
        a1 <- stats::aggregate( dat[,ii] , list( group) , mean , na.rm=T )[,2]
		h1$pR <- a1[1]
		h1$pF <- a1[2]
		h1$pdiff <- h1$pR - h1$pF
		h1$pdiff.adj <- NA	
		h1$uniformDIF <- mod2$coef[3]
		h1$se.uniformDIF <- sqrt( diag( stats::vcov(mod2)) )[3]
		h1$t.uniformDIF <- mod2$coef[3] / sqrt( diag( stats::vcov(mod2) ) )[3] 
		h1$sig.uniformDIF <- ""
		if ( h1$t.uniformDIF > quant ){ h1$sig.uniformDIF <- "+" }
		if ( h1$t.uniformDIF < - quant ){ h1$sig.uniformDIF <- "-" }
		h1$DIF.ETS <- ""
		#****
		# nonuniform DIF
		h1$nonuniformDIF <- mod3$coef[4]
		h1$se.nonuniformDIF <- sqrt( diag( stats::vcov(mod3)) )[4]
		h1$t.nonuniformDIF <- mod3$coef[4] / sqrt( diag( stats::vcov(mod3) ) )[4] 
		h1$sig.nonuniformDIF <- ""
		if ( h1$t.nonuniformDIF > quant ){ h1$sig.nonuniformDIF <- "+" }
		if ( h1$t.nonuniformDIF < - quant ){ h1$sig.nonuniformDIF <- "-" }
		matr <- rbind( matr , h1 )
        cat( ii , " " ) ; utils::flush.console()
        if ( ii %% 15 == 0 ){ cat("\n") }
        }
    cat("\n")
    # include variable of adjusted p values
#    ind <- which( colnames(matr) == "pF" )       
    matr[ , "pdiff.adj" ] <- matr$pR - matr$pF - mean( matr$pR - matr$pF  )   
	ind1 <- grep( "ETS" , colnames(matr) )
	#***
	# DIF ETS classifiaction
	stat <- abs( matr$uniformDIF )
	stat.low <- stat - quant * matr$se.uniformDIF
	matr[,"DIF.ETS"] <- "B"
	# DIF classification C
	ind <- which( ( stat > .64 ) & ( stat.low > .43 ) )
	if (length(ind) > 0){ matr[ind, "DIF.ETS"] <- "C" }
	# DIF classification A
	ind <- which( ( stat < .43 ) | ( stat.low < 0 ) )
	if (length(ind) > 0){ matr[ind, "DIF.ETS"] <- "A" }
	matr$DIF.ETS <- paste0( matr$DIF.ETS , 
			ifelse( matr$uniformDIF > 0 , "+" , "-" )	 )
	#*********************************
	# calculation of DIF variance
    dif1 <- dif.variance( dif=matr$uniformDIF , se.dif = matr$se.uniformDIF )		
	matr <- data.frame( matr[ , seq(1,ind1)] , "uniform.EBDIF" = dif1$eb.dif ,
		"DIF.SD" = dif1$unweighted.DIFSD , matr[ , seq(ind1+1 , ncol(matr)) ] )
	cat( paste0("\nDIF SD = " , round( 	dif1$unweighted.DIFSD , 3 ) ) , "\n")
	# sorting of the items
	g1 <- rank( matr$uniformDIF )
	matr <- data.frame( "itemnr" = 1:nrow(matr) ,
			"sortDIFindex" = g1 , 
				matr )
    # matr <- data.frame( "item" = colnames(dat) , matr )
    return(matr)
    }
#------------------------------------------------------------------------------

