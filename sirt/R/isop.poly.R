######################################################################
# Fitting the ISOP and ADISOP model
isop.poly <- function( dat , score.breaks=seq(0,1,len=10 ) , 
			conv=.0001 , maxit=1000 , epsilon=.025 ,
			progress=TRUE ){
	#*****************************************
	# convert matrix into matrix format
	dat <- as.matrix(dat)
	I <- ncol(dat)	
	# define response matrix
	dat.resp <- 1-is.na( dat )
    # scoring of persons
	res <- isop.scoring( dat )
	person <- res$person	
	item <- res$item
	p.itemcat <- res$p.itemcat
	score.itemcat <- res$score.itemcat	
	distr.fct <- res$distr.fct
	K <- max( dat , na.rm=TRUE)
	# person score is the modified percentile score
    stud.p <- person$mpsc
	# different scores of students
	if ( ! is.null( score.breaks) ){
	    qsc <- stats::quantile( stud.p , probs=score.breaks) 
		qsc[1] <- -1.001
		qsc[ length(qsc) ] <- 1.001
		qsc <- unique( qsc)
		stud.p <- cut( stud.p , breaks=qsc , labels=FALSE)
						}	
    person$scoregroup <- stud.p						
	# score groups students
	scores <- sort( unique( stud.p ) )
	SC <- length(scores)
	# scoring of students
    # create item groups
	# create table
	freq.correct <- array( NA , dim=c(SC , I , K) )
	freq.categories <- array( NA , dim=c(SC , I , K+1) )	
	dimnames(freq.categories)[[1]] <- 
			dimnames(freq.correct)[[1]] <- paste0( "score_" , scores )
	dimnames(freq.categories)[[2]] <-  
			dimnames(freq.correct)[[2]] <- colnames(dat)
	dimnames(freq.categories)[[3]] <- paste0( "Cat" , 0:K )			
	dimnames(freq.correct)[[3]] <- paste0( "Cat" , 0:(K-1) )
	# compute frequency table
	for (kk in 0:K){
		# kk <- 0
		if (kk < K ){
			a1 <- stats::aggregate( 1*(dat <= kk ) , list( stud.p ) , mean , na.rm=TRUE )
			freq.correct[,,kk+1] <- as.matrix(a1[,-1])		
					}
		a2 <- stats::aggregate( 1*(dat == kk ) , list( stud.p ) , sum , na.rm=TRUE )		
		freq.categories[,,kk+1] <- as.matrix(a2[,-1])		
			}
	
	# weights
	wgt <- freq.correct[,,1]
	a1 <- stats::aggregate( dat.resp , list( stud.p ) , sum , na.rm=TRUE )
	wgt <- a1[,-1]			
	freq.isop.fitted <- freq.correct.ordered <- list(1:K)
	for (kk in 0:(K-1) ){
		freq.correct.ordered[[kk+1]] <- 
			freq.correct[, order( item$pscore ) ,kk+1]    
							}
	######################################################						
	#***************
	# calculate ISOP Model (for all categories)
	for (kk in 0:(K-1) ){
		if (progress){ cat("\nFit ISOP Category" , kk  ) }
		# kk <- 1
		fc <- - freq.correct.ordered[[kk+1]]
		wc <- wgt[ , colnames(fc) ]
		colnames(wc) -> wcc
		res.isop <- fit.isop( freq.correct=fc , wgt=wc , conv=conv , 
				        maxit=maxit , progress=progress , calc.ll=FALSE)	
		freq.isop.fitted[[kk+1]] <- - res.isop$fX
					}
	##############################################################
	#*********	
	# ADISOP Model (start with the fitted isop solution)
	# compute survivor function
	surv.isop <- as.list( 1:K )
	for (kk in 0:(K-1)){ 
		surv.isop[[kk+1]] <- 1 - freq.isop.fitted[[kk+1]]
					}
	# compute sum of survivor functions
	sum.surv.isop <- surv.isop[[1]]
	for (kk in 2:K){
		sum.surv.isop <- sum.surv.isop + surv.isop[[kk]]
					}
	# ADISOP fit based on sum of survivor functions	
	res.adisop <- fit.adisop( freq.correct=sum.surv.isop / K , wgt , 
			conv=conv , maxit=maxit , 
			epsilon=epsilon , progress=progress  , calc.ll=FALSE )
	item.sc <- res.adisop$item.sc
	person.sc <- res.adisop$person.sc
	# second smoothing step in ADISOP
	# see Scheiblechner (2009)
	adfitted <- res.adisop$freq.fitted
	adfitted1 <- adfitted[,1:6]
	adfitted1 <- adfitted1[ order( 100000*adfitted1$item.index + adfitted1$stud.index ) , ]
	# create a long table for survivor function
	surv.smooth <- NULL

	for (kk in 1:K){
		surv.smooth.kk <- adfitted1
		surv.smooth.kk$categ <- kk
		surv.smooth.kk$survisop <- as.vector( surv.isop[[kk]] )
		surv.smooth <- rbind( surv.smooth , surv.smooth.kk )
					}
	surv.smooth$X1PX2 <- surv.smooth$item.sc + surv.smooth$person.sc
	surv.smooth <- surv.smooth[ order( 1000*rank( surv.smooth$X1PX2 )+
					100 - surv.smooth$categ )  , ]
	SV <- nrow(surv.smooth) / K 
	fc1 <- matrix( surv.smooth$survisop , nrow=K , ncol=SV )
	wc1 <- matrix( surv.smooth$wgt , nrow=K , ncol=SV )
	# fit isop
	res.isop2 <- fit.isop( freq.correct=fc1 , wgt=wc1 , calc.ll=FALSE ,
			progress=FALSE)
	r1 <- res.isop2$freq.fitted$freq.fitted
	surv.smooth$survadisop <- matrix( matrix( r1 , nrow=K , ncol=SV , byrow=TRUE ) , 
			ncol=1 )
	# compute distribution function
    surv.smooth$sm.distr.fct <- 1 - surv.smooth$survadisop	
	surv.smooth$item <- wcc[ surv.smooth$item.index ]
	F.orig <- array( NA , dim=c(I,K+1,SC ) )
	dimnames(F.orig)[[3]] <- paste0( "score_" , scores )
	dimnames(F.orig)[[1]] <- colnames(dat)
	dimnames(F.orig)[[2]] <- paste0( "Cat" , 0:K )	
	F.orig[,K+1,] <- 1
    F.adisop <- F.isop <- F.orig	
	F1 <- matrix( 1 , K+1 , SC )
	for (ii in 1:I){ 
		#	ii <- 1
		vv <- colnames(dat)[ii]
		surv.smooth.vv <- surv.smooth[ surv.smooth$item == vv , ]
		ind <- (surv.smooth.vv$item.index)[1]
		ind2 <-  as.matrix(surv.smooth.vv[ , c("categ","stud.index") ] )
		F1[ ind2 ] <- surv.smooth.vv$sm.distr.fct
		F.adisop[vv,,] <- F1
		F1[ ind2 ] <- 1-surv.smooth.vv$survisop
		F.isop[vv,,] <- F1		
				}
	F.sat1 <- aperm( freq.correct , c(2,3,1) )
	F.saturated <- F.orig
	F.saturated[ , 1:K , ] <- F.sat1
	##########################################
	#*******
	# fit graded response model
	res.grm <- fit.gradedresponse( freq.categories , SC , I , K , 
			conv=conv , maxit=maxit , progress=progress )
	
	#########################################################
	#*****
	# calculate probabilities
	prob.saturated <- prob.isop <- prob.adisop <- F.orig
	#***********
	# calculate probabilities
	prob.saturated <-  calc.prob.distr( Farray=F.saturated)
	prob.adisop <- calc.prob.distr( Farray=F.adisop)
	prob.isop <- calc.prob.distr( Farray=F.isop)
	prob.grm <- aperm( res.grm$prob , c(2,3,1) )
	# [ items , categories , score groups ]
	# calculate survivor functions P( X>=k)
	surv.saturated <- calc.survfct.isop(prob.saturated)
	surv.isop <- calc.survfct.isop(prob.isop)
	surv.adisop <- calc.survfct.isop(prob.adisop)
	surv.grm <- calc.survfct.isop(prob.grm)	
	
	#***
	# calculate log-likelihood
	# calculate frequencies
	freq.saturated <- prob.saturated
	for (kk in 1:(K+1) ){
		freq.saturated[ , kk , ] <- t(wgt) * prob.saturated[,kk,]
					}
	ll <- data.frame( "model" = c("saturated" , "isop" , "adisop" , "grm") , 
	   "ll" = rep(NA,4) )
	ll$ll[1] <- calc.ll.isop.poly( freq1=freq.saturated , prob1=prob.saturated )  	
	ll$ll[2] <- calc.ll.isop.poly( freq1=freq.saturated , prob1=prob.isop )	
	ll$ll[3] <- calc.ll.isop.poly( freq1=freq.saturated , prob1=prob.adisop )
	ll$ll[4] <- res.grm$ll
	NW <- mean( colSums(wgt) )	
	ll$llcase <- ll$ll / NW
	#*******
	#* collect results
	res <- list( "freq.correct"=freq.correct ,
			  "wgt"=wgt , 
			  "prob.saturated" = prob.saturated , 
			  "prob.isop" = prob.isop , "prob.adisop" = prob.adisop ,
			  "prob.grm"=prob.grm , 
			  "surv.saturated" = surv.saturated , 
			  "surv.isop" = surv.isop , "surv.adisop" = surv.adisop ,
			  "surv.grm"=surv.grm , 
			  "ll" = ll , 
			  "person"=person , "item"=item , 
			  "p.itemcat"=p.itemcat , "score.itemcat"=score.itemcat , 
			  "fit.grm" = res.grm , 
			  "s2" = Sys.time() , "dat" = dat ,
			  "model"="isop.poly")
	class(res) <- "isop"
	return(res)		
		}
##########################################################
# auxiliary function for conversion of a distribution
# function into a probability table
calc.prob.distr <- function(Farray){
		prob1 <- Farray
		KK <- dim(Farray)[2]
		prob1[,1,] <- Farray[,1,]
		for (kk in 2:KK){ prob1[,kk,] <- Farray[,kk,] - Farray[,kk-1,] }
        prob1[ prob1 < 0 ] <- 0
		prob1[ prob1 > 1 ] <- 1		
		return(prob1)
			}
############################################################
# calculate log-likelihood
calc.ll.isop.poly <- function( freq1 , prob1 , eps=10^(-14) ){
        prob1[ prob1 < 0 ] <- 0
		prob1[ prob1 > 1 ] <- 1
  		ll <- sum( freq1 * log( prob1 + eps ) , na.rm=TRUE)
		return(ll)
				}
#############################################################
# calculate survivor functions
calc.survfct.isop <- function(prob1){
		prob2 <- prob1
		K <- dim(prob1)[2]
		for ( kk in seq(K-1,1,-1) ){
			prob2[,kk,] <- prob2[ , kk+1 , ] + prob1[,kk,]
							}
		return(prob2)
			}