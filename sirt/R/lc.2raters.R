

##############################################################
# latent class model for two exchangeable raters
lc.2raters <- function( data , conv=.001 , maxiter=1000 , progress=TRUE ){
	#**********************
	data <- stats::na.omit(data)
	# data preparation and creation of a frequency table
	res <- lc2.data.prep(data)
	maxK <- res$maxK
	m1 <- res$m1
	
	m2.combs <- t( combinat::combn( maxK+1 , 2 ) )
	N2 <- nrow(m2.combs)
	m2.combs <- rbind( m2.combs , m2.combs[,c(2,1)] , cbind( 1:(maxK+1) , 1:(maxK+1 ) ) )
	m2.combs <- m2.combs[ order(paste( m2.combs[,1],m2.combs[,2])) , ]	
	m2 <- cbind( m2.combs , NA )
	m2[ , 3 ] <- m1[ m2[,1:2] ]
	rownames(m2) <- paste0( rownames(m1)[ m2[,1] ] , "-" ,  rownames(m1)[ m2[,2] ] )
	pi.k <- colSums( m1 )
	pi.k0 <- pi.k <- pi.k / sum( pi.k )
	TP <- maxK+1
	KK <- nrow(m2)
	# init probabilities
	probs <- matrix( 1/TP , nrow=TP , ncol=TP)
	rownames(probs) <- colnames(probs) <- rownames(m1)
	diag(probs) <- 1
	probs <- probs / rowSums( probs )
	counts <- probs
	par.change <- 1
	iter <- 0
	#********************
	#**** begin iterations
	while ( ( par.change > conv ) & ( iter < maxiter) ){
        probs0 <- probs	
		# compute likelihood
		f.yi.qk <- probs[ m2[,1] , ] *  probs[ m2[,2] , ]
		rownames(f.yi.qk) <- rownames(m2)
		# compute posterior
		f.qk.yi <- matrix( pi.k , nrow=KK , ncol=TP , byrow=TRUE ) * f.yi.qk
		f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )
		# compute expected counts
		nik <- f.qk.yi * m2[,3 ]
		Nk <- colSums( nik )
		pi.k <- Nk / sum( Nk )
		# compute updated probabilities
		for (zz in 1:TP){      # zz <- 2
			counts[zz,] <- colSums( m2[ , 3 ] * 1 * ( m2[,1] == zz ) * f.qk.yi )
					}
		# update probabilities
		probs <- counts / matrix( colSums( counts ) , TP , TP , byrow=TRUE )
		#---
		par.change <- max( abs(probs - probs0) )
		iter <- iter + 1
		if (progress){
		   cat( paste0("Iteration " , iter , " | Max. probability change = " ,
				round( par.change , 6 ) , "\n") )
			utils::flush.console()	
					}
			} # end iter
	#******* end iterations
	#*******************************										
	colnames(probs) <- paste0("Lat" , colnames(probs))	
	
	#*** posterior probabilities
	classprob.2raters.like <- f.yi.qk[ m2[,1] <= m2[,2] , ]
	classprob.2raters.like <- classprob.2raters.like / rowSums( classprob.2raters.like )

	classprob.2raters.post <- f.yi.qk[ m2[,1] <= m2[,2] , ]
	classprob.2raters.post <- classprob.2raters.post / rowSums( classprob.2raters.post )

	classprob.1rater.post <- classprob.1rater.like <- 0*probs
	
	for ( kk in 1:(maxK+1) ){
		f1 <- colSums( f.yi.qk[ ( m2[,1] == kk ) | ( m2[,2] == kk ) , , drop=FALSE] )
		classprob.1rater.like[kk,] <- f1
		classprob.1rater.like[kk,] <- classprob.1rater.like[kk,] / sum( classprob.1rater.like[kk,] )
		classprob.1rater.post[kk,] <- f1 * pi.k
		classprob.1rater.post[kk,] <- classprob.1rater.post[kk,] / sum( classprob.1rater.post[kk,] )				
					}	
					
	# arrange output
	f.yi.qk <- f.yi.qk[ m2[,1] <= m2[,2] , ]
	f.qk.yi <- f.qk.yi[ m2[,1] <= m2[,2] , ]
	
	m2 <- as.data.frame(m2)
	colnames(m2) <- c("cat1" , "cat2" , "AbsFreq" )
	m2$RelFreq <- m2$AbsFreq / sum( m2$AbsFreq )
	s2 <- Sys.time()
	
	# agreement statistics
	agree.stats <- lc2.agreement( m1 )
	
	#*** output list
	res <- list( "classprob.1rater.like" = classprob.1rater.like , 
	    "classprob.1rater.post" = classprob.1rater.post , 
		"classprob.2raters.like"= classprob.2raters.like , 
		"classprob.2raters.post"= classprob.2raters.post ,
		"f.yi.qk" = f.yi.qk , "f.qk.yi"=f.qk.yi , "probs" = probs ,
		"pi.k" = pi.k , "pi.k.obs" = pi.k0 , "freq.long" = m2 , "freq.table" = m1 ,
		"agree.stats" = agree.stats , 
		"data"=data , "N.categ" = maxK+1 , "s2"=s2)
	class(res) <- "lc.2raters"	
    return(res)	
	}
###############################################################
# 	S3 summary method for lc.2raters
summary.lc.2raters <- function( object , ... ){

	cat("---------------------------------------------------------------\n")
		d1 <- utils::packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
    cat("Latent Class Model for Two Exchangeable Raters \n")
	cat("---------------------------------------------------------------\n")
    cat( "Number of persons =" , nrow(object$data) , "\n" )    
    cat( "Number of item categories =" , object$N.categ , "\n" )    
	cat("\n********************************************\n")
	cat("Symmetrized Frequency Table\n")
	obji <- object$freq.table
	print(obji)
# include agreement statistics here	    

	cat("\n********************************************\n")
	cat("Rater agreement statistics\n\n")
    cat( "Percentage agreement =" , round(object$agree.stats["agree0"],4) , "\n" )    	
    cat( "Percentage agreement (up to 1 category) =" , round(object$agree.stats["agree1"],4) , "\n" )    	
    cat( "Cohen's kappa =" , round(object$agree.stats["kappa"],4) , "\n" )    		
    cat( "Weighted kappa (linear weights) =" , round(object$agree.stats["wtd.kappa.linear"],4) , "\n" )    			
    cat( "Gwet's AC1 =" , round(object$agree.stats["AC1"],4) , "\n" )    			
    cat( "Aickin's alpha =" , round(object$agree.stats["alpha.aickin"],4) , "\n" )    				
	cat("\n********************************************\n")
	cat("Item probabilities given latent classes\n")
	obji <- object$probs
	for (vv in 1:( ncol(obji))){
		obji[,vv] <- round( obji[,vv] , 4 )
					}
	print(obji)	
	cat("\n********************************************\n")
	cat("Latent class classification probability for one rater (likelihood estimate)\n")
	obji <- object$classprob.1rater.like
	for (vv in 1:( ncol(obji))){
		obji[,vv] <- round( obji[,vv] , 4 )
					}
	print(obji)		
	cat("\n********************************************\n")
	cat("Latent class classification probability for two raters (likelihood estimate)\n")
	obji <- object$classprob.2raters.like
	for (vv in 1:( ncol(obji))){
		obji[,vv] <- round( obji[,vv] , 4 )
					}
	print(obji)		
                }
#*******************************************************
