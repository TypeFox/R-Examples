#######################################
# attach all elements in an environment   
.attach.environment <- function( res , envir ){
#	e1 <- environment()
	CC <- length(res)
	for (cc in 1:CC){
		assign( names(res)[cc] , res[[cc]] , envir=envir )		
					}
			}

############################################
# gdm data preparation
.gdm.data.prep <- function( dat , data , weights , group ){
	I <- ncol(dat)
	N <- nrow(dat)
	data1 <- data
	data1[ is.na(data1) ] <- 9
	cat("************************************************************\n")
	cat("Data preparation\n")
	cat("Number of rows in data =" , nrow(data1) , "\n") ; utils::flush.console()
	item.patt.subj <- data1[,1]
	for ( ii in 2:I){
		item.patt.subj <- paste( item.patt.subj  , data1[,ii] , sep="")
				}
	#****
	# arrange groups
	if ( is.null(group)){ 
		G <- 1 
		group <- rep(1,N)
				} 
			else {
		gr2 <- unique( sort(paste( group ) ))
		G <- length(gr2)
		group <- match( group , gr2 )
							}	
    # calculate frequency of each item response pattern
	# case of one group
	if (G==1){ 
		if ( is.null(weights) ){ weights <- rep(1,N) }
		a2 <- rowsum( weights , item.patt.subj) 
		item.patt <- a2[,1]	
		# define data frame 'item.patt' with item response pattern and its frequency (weight)
		item.patt <- cbind( "pattern" = names(item.patt), 
				"freq" = as.numeric(as.vector( item.patt ) ) )
		weights <- as.numeric(paste(item.patt[,"freq"]))
			}	
	#***
	# multiple group case
	if ( is.null(weights) ){ 
			weights <- rep(1,N) 
				}
	if (G>1){
		for (gg in 1:G){
	#gg <- 1
		ind.gg <- which( group == gg )
		a2 <- rowsum( weights[ind.gg] , item.patt.subj[ind.gg] )
		a2 <- data.frame( "pattern" = rownames(a2) , a2[,1] )		
		colnames(a2)[2] <- paste0("freq.Gr" , gg)
		rownames(a2) <- NULL
		if (gg == 1){ item.patt <- a2 }
		if (gg > 1){
			item.patt <- merge( item.patt , a2 , by ="pattern" , all=TRUE )
					}
		item.patt[ is.na(item.patt) ] <- 0
				}
		weights <- item.patt[,-1]
				}	
# print(item.patt)	
	#***
	# reconstruct data
	N <- nrow(item.patt)
	dat <- matrix(NA , N , I )
	for (ii in 1:I){
		dat[,ii ] <- as.numeric( substring( item.patt[,"pattern"] , ii,ii) )
				}
	dat.resp <- 1-(dat==9)
	data <- dat
	dat[ dat.resp==0] <- 0
	cat("Number of response patterns =" , nrow(dat) , "\n")
	utils::flush.console()
	res <- list( "weights"=weights , "dat" = dat , "dat.resp"=dat.resp ,
		"data"=data , "item.patt" = item.patt )
	return(res)
	}			
####################################################		
			
#################################################
# define Q matrix
.gdm.Qmatrix <- function(Qmatrix,irtmodel,I,TD,K,a){
	# Q matrix [1:I , 1:TD , 1:K]
	if ( is.null(Qmatrix) ){
		Qmatrix <- array( 1 , dim=c(I,TD,K) )	
			# modify it possibly
		if (K>1 & ( irtmodel != "2PLcat" ) ){
			for (kk in 2:K){Qmatrix[,,kk] <- kk*Qmatrix[,,1] }
				}				
		if ( irtmodel=="2PLcat"){
			for (kk in 2:(K) ){
				a[,,kk] <- kk * a[,,kk]
				}
			}
		}
	res <- list("Qmatrix" = Qmatrix , "a" = a)
	return(res)
		}			
			
##########################################
# Theta design matrix
.gdm.thetadesign <- function( theta.k , thetaDes , Qmatrix ){
	
	D <- 1  # default dimension 1
	####################
	# definition of theta.k
	if ( ! is.null(Qmatrix) ){
		D <- ncol(Qmatrix)
			if ( length( dim(Qmatrix))==2 ){ 
				Q1 <- array( 0 , dim=c(dim(Qmatrix),1) )
				Q1[,,1] <- Qmatrix
				Qmatrix <- Q1
						}
						}
	w1 <- ( is.vector( theta.k)  ) & ( ! is.list( theta.k) )	
	if ( w1 ){
#	if ( is.numeric(theta.k) ){
			theta.k <- matrix( theta.k , ncol=1 )
			if (D>1){
				th1 <- as.list(1:D)
				for (dd in 1:D){ th1[[dd]] <- theta.k }
				theta.k <- th1
					}
						}	
	if ( is.list( theta.k) ){
			tk <- theta.k	
			theta.k <- expand.grid( theta.k )
			colnames(theta.k) <- names(tk)
								}
							
	theta.k <- as.matrix(theta.k)	
	D <- ncol(theta.k)											
	if ( is.null( colnames(theta.k) ) ){
		colnames(theta.k) <- paste0("F",1:D)
				}		
	##############################																							
	if ( is.null(thetaDes) ){
		# thetaDes [TP,TD]
		TD <- D
		thetaDes <- matrix( theta.k , ncol=TD )
		colnames(theta.k) -> colnames(thetaDes)
					}
#	if (is.null( colnames(thetaDes) ) ){
#		colnames(thetaDes) <- paste0( "F" , 1:TD )
#					}
	TP <- nrow(thetaDes)	
	TD <- ncol(thetaDes)
	res <- list("D"=D , "TD"=TD , "TP"=TP,"theta.k"=theta.k,
		"thetaDes"=thetaDes , "Qmatrix"=Qmatrix )
	return(res)
		}

########################################################
# create delta design matrix
.gdm.create.delta.designmatrix <- function( delta.designmatrix , 
		TP , D , theta.k , skill.levels,G){
	delta.designmatrix <- rep(1,TP)
	for (dd in 1:D){
		# dd <- 1
		for ( pp in 1:(min( skill.levels[dd]-1 ,3) ) ){
			delta.designmatrix <- cbind( delta.designmatrix , theta.k[,dd]^pp )
						}
					}
	if (D>1){
		for (dd1 in 1:(D-1) ){				
			for (dd2 in (dd1+1):D) {					
				delta.designmatrix <- cbind( delta.designmatrix , theta.k[,dd1]*theta.k[,dd2] )		
								}
							}		
						}
	delta <- matrix(0,ncol(delta.designmatrix),G)
	covdelta <- NULL
			
	res <- list( "delta" = delta , "covdelta" = covdelta , 
		"delta.designmatrix" = delta.designmatrix )
	return(res)
		}
		
		
###############################################
# constraints for item parameters
.gdm.constraints.itempars <- function( b.constraint , a.constraint , 
	K , TD , Qmatrix , a ){
	for (kk in 1:K){
	  for( td in 1:TD){
		# kk <- 1 ; td <- 1
			ind.kk <- which( Qmatrix[ ,td , kk] == 0 )
			a[ ind.kk , td , kk ] <- 0
			if ( length( ind.kk) > 0 ){
				a1 <- cbind( ind.kk  , td , kk , 0 )
				a.constraint <- rbind( a.constraint , a1 )
							}
						}
					}
	if ( ! is.null( a.constraint) ){
		a.constraint <- as.matrix( a.constraint )
					}				
    res <- list( "a.constraint" = a.constraint , "b.constraint"=b.constraint , "a"=a)			
	return(res)
		}

##########################################################
# constraints on item parameters
.gdm.constraints.itempars2 <- function( b.constraint , a.constraint , 
	K , TD ,I , dat ){
	K.item <- apply( dat , 2 , max )
	for (ii in 1:I){	# ii <- 1
	K.ii <- K.item[ii]	
	if ( K.ii < K ){
		for ( kk in (K.ii+1):K){
			b.constraint <- rbind( b.constraint , cbind( ii , kk , -99999 ) )
		   for (td in 1:TD){
#				a.constraint <- rbind( a.constraint , cbind( ii , td , kk , -99999 ) )
				a.constraint <- rbind( a.constraint , cbind( ii , td , kk , 0 ) )
							}
							}
						}
				}
	res <- list("K.item"=K.item , "a.constraint"=a.constraint ,
			"b.constraint" = b.constraint )
	return(res)
	}
###############################################################			