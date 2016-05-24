
###################################################
# data preparation function
# linear logistic partial credit model
lpcm_data_prep <- function( dat , weights , a ){
    used_persons <- as.vector(which(rowSums( 1 - is.na(dat) ) > 1 ))
	dat <- dat[ used_persons , ]
	if ( ! is.null(weights) ){
		weights <- weights[ used_persons ]
							}
	N <- nrow(dat)
	I <- ncol(dat)
	if ( is.null(a) ){
		a <- rep(1,I)
					}
	aM <- matrix( a , nrow=N , ncol=I , byrow=TRUE )
	dat <- dat * aM
	
	dat_ind <- 1*is.na(dat)
	res <- sirt::md.pattern.sirt(dat_ind)
	resp_patt <- res$resp_patt
	
	
    if ( is.null(weights) ){
		weights <- rep(1 , N )
							}
	patt_unique <- unique( resp_patt )
	patt <- match( resp_patt , patt_unique )
	NP <- length(patt_unique)
	maxK <- apply( dat , 2 , max , na.rm=TRUE ) 
	item_index <- sapply( 1:NP , FUN = function(pp){
			# pp <- 1
			dat0 <- dat_ind[  patt == pp , , drop=FALSE ][1,]
			which( dat0 == 0 )
					} , simplify = FALSE )
	maxscore <- sapply( 1:NP , FUN = function(pp){
		sum(maxK[ item_index[[pp]] ] )
			} , simplify=FALSE )
	pars <- rep( 1:I , maxK )		
	pars_info <- data.frame( "item" = rep(colnames(dat),maxK) , "itemid" = pars )
	pars_info$cat <- unlist( sapply( 1:I, FUN = function(ii){
				seq( 1 , maxK[ii] ) 
						}  , simplify=FALSE) )
	pars_info$index <- seq(1 , nrow(pars_info) )
	pars_info$maxK <- maxK[ pars_info$itemid ]
	pars_info$Freq <- 0
	
	K <- max( maxK )
	for (kk in 1:K){
    	f1 <- colSums( ( dat == kk ) * weights , na.rm=TRUE)
		ind <- pars_info[ pars_info$cat == kk , "itemid" ]
		pars_info[ pars_info$cat == kk , "Freq" ] <- f1[ind]
					}
	parm_index <- sapply( 1:NP , FUN = function(pp){
		 which( pars_info$itemid %in% item_index[[pp]] ) 
			} , simplify=FALSE )
			
			
	# calculate raw score
	rs <- rowSums( dat , na.rm=TRUE)
	score_freq <- sapply( 1:NP , FUN = function( pp ){
			# pp <- 2
			sapply( 0:maxscore[[pp]] , FUN = function(ss){
								sum( ( rs == ss ) * ( patt == pp ) * weights , na.rm=TRUE) } ,
								simplify = TRUE )
						} , simplify=FALSE )

		suffstat <- sapply( 1:NP , FUN = function(pp){
			unlist( sapply( item_index[[pp]] , FUN = function(ii){
						sapply( 1:maxK[ii] , FUN = function(kk){
								sum( ( dat[,ii] == kk ) * weights * ( patt == pp ) , na.rm=TRUE) 
												} , simplify=FALSE )
											}
							) ) } , simplify = FALSE )

		splitvec <- sapply( 1:NP , FUN = function(pp){
				mv <- item_index[[pp]] 
				oj_max <- maxK[ mv ]
				rep.int( mv , oj_max)  
						} , simplify=FALSE )

	# generate pararameter names
	parnames <- paste0( pars_info$item , "_Cat" , pars_info$cat )	
	rownames(pars_info) <- parnames					
						
	res <- list( N=N , I=I , NP=NP , dat=dat , 
			      patt=patt , weights=weights ,
				  suffstat=suffstat, splitvec=splitvec,
				  item_index = item_index , parm_index=parm_index ,
				  pars_info=pars_info , maxscore=maxscore ,
				  maxK=maxK , score=rs , used_persons = used_persons ,
				  score_freq = score_freq , a = a , parnames = parnames )	
	return(res)	
		    
		}
###################################################