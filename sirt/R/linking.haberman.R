
##################################################
# Linking Haberman ETS Research Report
linking.haberman <- function( itempars , personpars=NULL ,
	conv = .00001 , maxiter=1000 ,progress=TRUE ){
	#****
    # convert itempars to data frame
	itempars <- as.data.frame( itempars )
	# include wgt if there does not exist a fifth columm
	if ( ncol(itempars) == 4){
		itempars$wgt <- 1
				}
	# extract studies
	studies <- sort( paste( unique( itempars[,1] ) ) )
	NS <- length(studies)
	# extract items
	items <- sort( paste( unique( itempars[,2] ) ) )
	NI <- length(items)
	# define a and b matrices
	wgtM <- bM <- aM <- matrix(NA , nrow=NI , ncol=NS)
	rownames(wgtM) <- rownames(bM) <- rownames(aM) <- items
	colnames(wgtM) <- colnames(bM) <- colnames(aM) <- studies
	# define item parameters
	for (ss in studies){
		# ss <- studies[1]
		itempars.ss <- itempars[ itempars[,1] == ss , ]
		aM[ paste(itempars.ss[,2]) , ss ] <- itempars.ss[,3] 
		bM[ paste(itempars.ss[,2]) , ss ] <- itempars.ss[,4] 
		wgtM[ paste(itempars.ss[,2]) , ss ] <- itempars.ss[,5] 
					}
	a.orig <- aM
	b.orig <- bM
	wgtM <- wgtM / matrix( rowSums( wgtM , na.rm=TRUE ) , nrow=NI , ncol=NS )
	#*****
	# estimation of A
	logaM <- log( aM )	
	resA <- .linking.haberman.als(logaM=logaM , wgtM=wgtM , maxiter=maxiter , 
				 conv=conv , progress = progress , est.type="A (slopes)")
	aj <- exp( resA$logaj )
	At <- exp( resA$logaAt )
	#******
	# estimation of B
	#	bMadj <- bM / matrix( At , NI , NS , byrow=TRUE )
	# correction ARb 2013-10-09
	bMadj <- bM * matrix( At , NI , NS , byrow=TRUE )
	resB <- .linking.haberman.als(logaM=bMadj , wgtM=wgtM , maxiter=maxiter , 
				 conv=conv , progress = progress , est.type="B (intercepts)")
	Bj <- resB$logaj
	Bt <- resB$logaAt
	#*****
	# transformations
	transf.pars <- data.frame( "study" = studies , "At" = At , "Bt" = Bt )
	rownames(transf.pars) <- NULL	
	transf.itempars <- data.frame( "study" = studies , "At" = 1/At , 
			"At2"=At , "Bt" = Bt )
	rownames(transf.itempars) <- NULL
	# This is the transformation for item parameters.
	#****
	transf.personpars <- transf.itempars[,c(1,2,4)]
	transf.personpars$At <- transf.pars$At
#	transf.personpars$Bt <-  - transf.itempars$Bt / transf.itempars$At	
	transf.personpars$Bt <-  - transf.pars$Bt 
#	transf.personpars <- transf.personpars
	colnames(transf.personpars) <- c("study" , "A_theta" , "B_theta" )
	colnames(transf.itempars) <- c("study" , "A_a" , "A_b" ,"B_b" )
	# new item parameters
	joint.itempars <- data.frame("item" = items , "aj"=aj , "bj"=Bj )
	# transformations for item parameters
	aM <- aM / matrix( At , NI , NS , byrow=TRUE )
	bM <- bM * matrix( At , NI , NS , byrow=TRUE ) - matrix( Bt , NI , NS , byrow=TRUE )	
	#****
	# transform person parameters
	if ( ! is.null( personpars) ){
	  for (ll in 1:NS){
		pp0 <- pp1 <- personpars[[ll]]
		pp1 <- transf.personpars$A_theta[ll] * pp1 + transf.personpars$B_theta[ll]
		ind <- which( substring( colnames(pp0),1,2) %in% c("se" , "SE") )
		if ( length(ind) > 0 ){
				pp1[,ind] <- transf.personpars$A_theta[ll] * pp0[,ind]
							}
		ind <- which( substring( colnames(pp0),1,3) %in% c("pid") )
		if ( length(ind) > 0 ){
				pp1[,ind] <- pp0[,ind]
							}
		personpars[[ll]] <- pp1
					}
	  }
	#*************************
	# calculate R-squared measures of invariance

	# select items for R2 calculation for which at least
	# two studies are available.
	selitems <- which( rowSums( 1 - is.na( a.orig ) ) > 1 )
	
	Rsquared.invariance <- c(NA,NA)
	names(Rsquared.invariance) <- c("slopes" , "intercepts" )	
	
	# retransformed parameters
	aj1 <- aj * matrix( At , NI , NS , byrow=TRUE )
	a.res <- a.orig - aj1
	# check correctness of formulas!
	
	Rsquared.invariance["slopes"] <- 1 -
		sum( a.res[ selitems,]^2 , na.rm=TRUE ) / 
		sum( a.orig[ selitems , ]^2 , na.rm=TRUE )
#	Rsquared.invariance["slopes"] <- 1 -
#		sum( a.res[ selitems,]^2 , na.rm=TRUE ) / 
#		var( as.vector(a.orig[ selitems , ]^2) , na.rm=TRUE )
	bj1 <- 1 / matrix( At , NI , NS , byrow=TRUE )*( 
			Bj + matrix( Bt , NI , NS , byrow=TRUE ) )
	b.res <- b.orig - bj1
	Rsquared.invariance["intercepts"] <- 1 -
		sum( b.res[ selitems,]^2 , na.rm=TRUE ) / 
		sum( b.orig[ selitems , ]^2 , na.rm=TRUE )
#	Rsquared.invariance["intercepts"] <- 1 -
#		sum( b.res[ selitems,]^2 , na.rm=TRUE ) / 
#		var( as.vector(b.orig[ selitems , ]^2) , na.rm=TRUE )
	es.invariance <- rbind( Rsquared.invariance ,
			sqrt( 1- Rsquared.invariance ) )
	rownames(es.invariance) <- c("R2" , "sqrtU2")
	
    # print output
	if (progress){
	    cat("Transformation parameters (Haberman linking)\n")
		obji <- transf.pars
		obji[,-1] <- round( obji[,-1] , 3)
		print( obji ) 	
	    cat("\nLinear transformation for item parameters a and b\n")
		obji <- transf.itempars
		obji[,-1] <- round( obji[,-1] , 3)
		print( obji )  	
	    cat("\nLinear transformation for person parameters theta\n")
		obji <- transf.personpars
		obji[,-1] <- round( obji[,-1] , 3)
		print( obji )  			
	    cat("\nR-Squared Measures of Invariance\n")
		obji <- es.invariance
		obji <- round( obji , 4)
		print( obji )  			
				}
	res <- list( 
		"transf.pars" = transf.pars , 
		"transf.itempars" = transf.itempars , 
	    "transf.personpars" = transf.personpars , 
		"joint.itempars" = joint.itempars ,
		"a.trans" = aM ,
		"b.trans" = bM ,
		"a.orig" = a.orig , "b.orig" = b.orig , 
		"personpars"=personpars,
		"es.invariance"=es.invariance)
	return(res)
		}
#######################################################################		
	
	