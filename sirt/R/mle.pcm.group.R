#################################################################
# MLE of person parameters for individuals and groups
mle.pcm.group <- function( dat , b , a =rep(1,ncol(dat)) ,
	group=NULL , pid = NULL , adj_eps = .3 , conv=.0001 ,
    maxiter=30	){
	# indicator for groupwise estimation
	est_group <- ! is.null( group )
	if ( is.null( pid )){ pid <- 1:(nrow(dat)) }
	#******
	# data preparation
	data <- dat
	dat_resp <- 1 - is.na(data)
	dat <- data
	dat[ is.na(data) ] <- 0
	dat <- as.matrix(dat)
	theta0 <- rep( 0 , nrow(dat) )
	dat_resp <- as.matrix( dat_resp )
	#*****
	# group indicators
	if (est_group){
#		group <- paste( group )
		ind <- order(group)
		dat_resp <- dat_resp[ ind , ]
		dat <- dat[ind,]
		group <- group[ind]
		idgroup <- match( group , unique( group ))
		# group indices
		gr1 <- rowsum( 1+0*idgroup , idgroup )
		c1 <- cumsum( gr1[, 1] ) 
		groupM <- as.matrix( cbind( c( 1 , c1[ - length(c1) ] +1 ) , c1 ) )
					}
	#*********************
	# extract maximum values
	maxK <- apply( dat , 2 , max )
	# compute person scores
	personScore <- rowSums( dat * dat_resp ) 
	# compute person maximum score
	personMax <- rowSums( matrix( maxK , nrow=nrow(dat) , ncol=ncol(dat) , byrow=TRUE ) * dat_resp )

	#******
	# epsilon adjustment for groups
	if ( est_group){
		groupMax <- rowsum( personMax , idgroup )
		groupScore <- rowsum( personScore , idgroup )
		grstat <- data.frame( groupScore , groupMax )
		fac <- ( groupMax - 2*adj_eps ) / groupMax
		maxKM <- dat_resp*matrix( maxK , nrow(dat) , ncol(dat) , byrow=TRUE )
		# adjustment matrix
		adjM <- maxKM / groupMax[ idgroup ] * adj_eps * dat_resp
		dat2 <- adjM + dat_resp * ( dat * fac[ idgroup ] )
		grstat$groupScoreAdj <- rowsum( rowSums( dat2*dat_resp ) , idgroup )
		theta0 <- rep( 0 , nrow(groupM) )
				}
	if ( ! est_group ){ # person estimation 
		fac <- ( personMax - 2*adj_eps ) / personMax
		maxKM <- dat_resp*matrix( maxK , nrow(dat) , ncol(dat) , byrow=TRUE )
		# adjustment matrix
		adjM <- maxKM / personMax * adj_eps * dat_resp
		dat2 <- adjM + dat_resp * ( dat * fac )
		personScoreAdj <- rowSums( dat2 * dat_resp ) 
	    theta0 <- rep(0, nrow(dat2) ) 
			}
	b <- as.matrix(b)
	dat <- dat2

	if ( est_group){
		res1 <- mle_pcm_group_CR( dat , dat_resp , groupM , b , a , 
				maxK+1 , theta0 , conv , maxiter )
		res <- data.frame( "group" = unique(group) , "idgroup" = rownames(gr1) ,
				"N.group" = gr1[,1] , 
				grstat , "theta" = res1$theta , "setheta"=res1$setheta ,
				"Niter" = res1$Niter )
				}
	if ( ! est_group){
		res1 <- mle_pcm_CR( dat , dat_resp , b , a , 
				maxK+1 , theta0 , conv , maxiter )
		res <- data.frame( "pid" = pid , "personScore"=personScore ,
				"personMax"=personMax , "personScoreAdj"=personScoreAdj ,
				"theta" = res1$theta , "setheta"=res1$setheta ,
				"Niter" = res1$Niter )
				}
	#*********
	# output
	res <- list( "person" = res , "data_adjeps" = dat )
	return(res)
}
########################################################
# Partial Credit Model: groupwise likelihoods
# extern "C" {
# SEXP mle_pcm_group_C( SEXP dat_, SEXP dat_resp_, SEXP groupM_, 
# SEXP b_, SEXP a_, SEXP maxK_, SEXP theta0_, SEXP conv_, SEXP maxiter_) ;
# }
mle_pcm_group_CR <- function (dat_, dat_resp_, groupM_ , b_, 
       a_, maxK_, theta0_, conv_, maxiter_ ){ 
	.Call("mle_pcm_group_C", dat_, dat_resp_, groupM_ , b_, 
			a_, maxK_, theta0_, conv_, maxiter_ , 
			PACKAGE = "sirt")
					}
########################################################


########################################################
# Partial Credit Model: Individual likelihoods
# extern "C" {
# SEXP mle_pcm_C( SEXP dat_, SEXP dat_resp_, SEXP b_, 
# SEXP a_, SEXP maxK_, SEXP theta0_, SEXP conv_, SEXP maxiter_) ;
# }
mle_pcm_CR <- function(dat_, dat_resp_, b_, 
       a_, maxK_, theta0_, conv_, maxiter_ ){ 
	.Call("mle_pcm_C", dat_, dat_resp_, b_, 
			a_, maxK_, theta0_, conv_, maxiter_ , 
			PACKAGE = "sirt")
					}
########################################################