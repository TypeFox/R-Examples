
##########################################################
# vcov din object
vcov.din <- function( object , extended = FALSE , infomat=FALSE ,
			ind.item.skillprobs = TRUE , ind.item= FALSE , diagcov = FALSE , h=.001 , ... ){

	#********			
# a0 <- Sys.time()			
    infomat.ind <- infomat
	latresp <- object$control$latresp
	guess <- object$guess$est
	slip <- object$slip$est

	item.patt.split <- object$item.patt.split
	item.patt.freq <- object$item.patt.freq
	attribute.patt <- object$attribute.patt$class.prob

	# calculate log-likelihood for every case
	weights <- item.patt.freq
	resp.ind.list <- object$control$resp.ind.list
	partable <- object$partable
	NPars <- max(partable$parindex)
	IP <- length(item.patt.freq)
	parnames <- unique( partable$parnames[ partable$parindex > 0 ] )
	ll.derivM <- matrix( NA , nrow=IP , ncol=NPars )
	colnames(ll.derivM) <- parnames
    
	J <- length(guess)
    L <- length(attribute.patt) 	

	#*** LL evaluated at theta
	guess0 <- guess
	slip0 <- slip
	skillprobs0 <- attribute.patt
	res <- vcov.loglike.din( weights , skillprobs0 , slip0 , guess0 ,
					latresp , item.patt.split , resp.ind.list , return.p.xi.aj=TRUE)
	ll1 <- res$ll
	p.xi.aj <- res$p.xi.aj
	
	#-----------------------------------------
	# compute first derivative
	for (pp in 1:NPars){
	
		#*** LL evaluated at theta
		guess0 <- guess
		slip0 <- slip
		skillprobs0 <- attribute.patt
		
		#*** LL evaluated at theta + h
		partable.pp <- partable[partable$parindex == pp ,]
		guess0 <- guess
		slip0 <- slip
		skillprobs0 <- attribute.patt
		recalc.ll <- TRUE
		if ( paste(partable.pp$partype[1]) == "guess" ){
			ind <- partable.pp$varyindex
			vecadd <- rep(0,J)
			vecadd[ind] <- h
			guess0 <- guess0 + vecadd 			
						}
		if ( paste(partable.pp$partype[1]) == "slip" ){
			ind <- partable.pp$varyindex
			vecadd <- rep(0,J)
			vecadd[ind] <- h
			slip0 <- slip0 + vecadd 
						}
		if ( paste(partable.pp$partype[1]) == "probs" ){
			ind <- partable.pp$varyindex
			vecadd <- rep(0,L)
			vecadd[ind] <- h
			skillprobs0 <- skillprobs0 + vecadd 
			recalc.ll <- FALSE            
						}                
		if ( recalc.ll){
			ll2 <- vcov.loglike.din( weights , skillprobs0 , slip0 , guess0 ,
						latresp , item.patt.split , resp.ind.list)
						} else {
			skillprobsM <- matrix( skillprobs0 , nrow=IP , ncol=L , byrow=TRUE )
			ll2 <- log( rowSums( p.xi.aj * skillprobsM ) )
						}
					
		ll.deriv1 <- ( ll2 - ll1 ) / h 
		ll.derivM[,pp] <- ll.deriv1	
		}
		
#cat("compute first derivative") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
		
	#----------------------------------------------------------------
	# compute information matrix
#	infomat <- matrix( 0 , nrow=NPars , ncol=NPars )
#	rownames(infomat) <- colnames(infomat) <- parnames

	
# cat("compute second derivative") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		

	#*************************
	infomat2 <- crossprod( ll.derivM , weights * ll.derivM )
	if (ind.item.skillprobs){
		pt1 <- partable[ partable$partype %in% c("guess","slip") , ]
		pt2 <- partable[ partable$partype %in% c("probs") , ]	
		pt1 <- unique(setdiff( pt1$parindex , 0 ))
		pt2 <- unique(setdiff( pt2$parindex , 0 ))
		infomat2[ pt1 , pt2 ] <- infomat2[ pt2 , pt1 ] <- 0
							}
	if (ind.item){
		pt1 <- partable[ partable$partype %in% c("guess","slip") , ]
		h1 <- expand.grid( pt1$parindex , pt1$parindex )
		h1 <- merge( x = h1 , y = pt1[ , c("parindex" , "item") ] , by = 1 )
		h1 <- merge( x = h1 , y = pt1[ , c("parindex" , "item") ] , by.x = 2 , by.y=1)
		h1 <- h1[ h1[,3] != h1[,4] , ]
		infomat2[ as.matrix( h1[,1:2] )  ] <- 0
				}							
						
# cat("compute second derivative (II)") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							
	infomat <- infomat2
	# inverse of information matrix
    if (infomat.ind ){
		covmat <- infomat 
		cat("The information matrix is computed.\n")
				} else {
		covmat <- solve( infomat )
		attr(covmat , "coef") <- coef(object)
		# extended set of coefficients
		if (extended){
			A <- object$vcov.derived$A
			covmat <- A %*% covmat %*% t(A)
			x1 <- object$partable$value
			names(x1) <- object$partable$parnames
			attr(covmat,"coef") <- x1
					}
				}						
	return(covmat)	
		}
#########################################################################
