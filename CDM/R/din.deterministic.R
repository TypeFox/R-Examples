
#####################################################
# Deterministic din estimation
din.deterministic <- function( dat , q.matrix , rule="DINA" , method="JML" , 
	conv=.001 , maxiter=300 , increment.factor=1.05 , progress=TRUE){
	#********************************
	# data preparations
	dat0 <- dat
	I <- ncol(dat)
	N <- nrow(dat)
	dat[ is.na(dat) ] <- 0
	dat.resp <- 1*(1-is.na(dat0))
	if ( length(rule) == 1 ){ rule <- rep( rule , I ) }
	dat <- as.matrix(dat)
	dat.resp <- as.matrix(dat.resp)
	# define attribute patterns
	attr.patt <- define.attribute.space( q.matrix )
	AP <- nrow(attr.patt)
	rownames(q.matrix) <- colnames(dat)
	# compute latent responses
	latresp <- compute.latent.response( attr.patt , q.matrix , rule)
	# initial guessing and slipping parameters
	guess <- stats::runif(  I  , .1 , .15)
	slip <- stats::runif(  I  , .1 , .15)
	max.increment <- 1
	# initial latent response vector
	latresp.est <- latresp[ rep(1,N) , ]
    parchange <- 1000
	if ( method=="weighted.hamming"){
		pbar <- colMeans( dat0 , na.rm=TRUE )
		guess <- slip <- 1 / ( pbar * ( 1 - pbar ) )
		maxiter <- 1
				}
	if ( method=="hamming"){
		guess <- slip <- rep(.2,I)
		maxiter <- 1
				}				
    iter <- 0
	while( ( iter < maxiter ) & (parchange > conv) ){
		slip0 <- slip
		guess0 <- guess
		latresp.est0 <- latresp.est	   
		# compute individual deviations
		if (method!="JML"){
			res <- din.deterministic.devcrit( dat=dat , datresp=dat.resp , latresp , guess , slip )
						  }
	    if (method=="JML"){
			res <- din.jml.devcrit( dat=dat , datresp=dat.resp , latresp , guess , slip )
						}
		# compute individual classifications
		attr.est <- attr.patt[ res$indexcrit , ]
		attr.est.index <- res$indexcrit
		# extract estimated latent responses
		latresp.est <- latresp[ attr.est.index , ]	
		# changes in estimated skill patterns
		latresp.change <- mean( abs( latresp.est - latresp.est0 ) )
		# calculate guessing and slipping parameters
		# slipping
		L1 <- ( latresp.est == 1 ) * (dat.resp==1 )
		L2 <- L1*(dat==1)
		slip <- 1 - colSums(L2) / ( colSums(L1) + .00001 )
		# guessing
		L1 <- ( latresp.est == 0 ) * (dat.resp==1 )
		L2 <- L1*(dat==0)
		guess <- 1 - colSums(L2) / ( colSums(L1) + .00001 )
		# update increment	
		if (maxiter > 1){
			increment <- guess - guess0
			max.increment <- max.increment/increment.factor
			guess <- guess0 +  ifelse( abs( increment) > max.increment , sign(increment)*max.increment , increment )
			increment <- slip - slip0
			slip <- slip0 +  ifelse( abs( increment) > max.increment  , sign(increment)*max.increment , increment )		
						}
			max.increment <- parchange <- max( abs( c(guess-guess0 , slip-slip0) ) )
		# extract deviation value
		if ( method!="JML"){
			devval <- sum(res$mincrit) } else {
			devval <- sum(-2*log(res$mincrit))	}
		iter <- iter + 1
		# print progress
		if (progress){
		    cat("Iteration" , iter )
			if ( method=="JML"){
				cat(" | Deviance = ",round(devval,3) )
				cat("\n ****")
								}
			if ( method!="JML"){ cat(" |") }								
				cat("  Average change in classifications = ",round(latresp.change,5) )
			if (method %in% c("JML","adaptive") ){	
					cat(" | Max. param. change = " , round( max.increment,6) , "\n")
												} else { cat("\n") }
			utils::flush.console()			
					}				
			}
		##################################################################
		# compute prediction error
		prederror <- sum( dat.resp * abs( dat - latresp.est ) ) / sum( dat.resp )
		if (progress){
		    cat("-------------------\n" )
			cat("Average predictor error=",round(prederror,9) )
				cat("\n")
								}		
		# collect output values			
		res <- list( "attr.est" = latresp.est , "criterion" = devval ,
			"guess"=guess , "slip"=slip , "prederror" = prederror , 
			"q.matrix" = q.matrix ,
			"dat" = dat0 )
	return(res)		
		}
#################################################################################