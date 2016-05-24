######################################################################
# Fitting the ISOP and ADISOP model
isop.dich <- function( dat , score.breaks=NULL , merge.extreme=TRUE ,
			conv=.0001 , maxit=1000 , epsilon=.025 ,
			progress=TRUE ){
	#*****************************************
	# convert matrix into matrix format
	dat <- as.matrix(dat)
	I <- ncol(dat)	
	# define response matrix
	dat.resp <- 1-is.na( dat )
	if( max( dat , na.rm=TRUE ) > 1 ){
		stop("Use isop.poly for polytomous data!")
					}
	if ( ! is.null(score.breaks) ){
		merge.extreme <- FALSE }
	# different scores of students
	stud.p <- rowMeans( dat , na.rm=TRUE )
	if (merge.extreme){
		score.breaks <- ( c(0,2:(I-1) , I+1 ) -.5 ) / I
					  }	
	if ( ! is.null( score.breaks) ){
		stud.p <- cut( stud.p , breaks=score.breaks , labels=FALSE)
						}	
	# different items
	item.p <- colMeans( dat , na.rm=TRUE )
	item.ps <- sort( item.p, index.return=TRUE)
	dat <- dat[ ,  item.ps$ix ]
	# score groups students
	scores <- sort( unique( stud.p ) )
	SC <- length(scores)
	# scoring of students
	res <- isop.scoring( dat )
	person <- res$person
	person$scoregroup <- stud.p
	item <- res$item
	p.itemcat <- res$p.itemcat
	score.itemcat <- res$score.itemcat
    # create item groups
	# create table
	freq.correct <- matrix( NA , SC , I )
	rownames(freq.correct) <- paste0("score_" , scores )
	colnames(freq.correct) <- colnames(dat)
	wgt <- freq.correct
	# percent correct
	a1 <- stats::aggregate( dat == 1 , list( stud.p ) , mean , na.rm=TRUE )
	freq.correct <- a1[,-1]
	# weights
	a1 <- stats::aggregate( dat.resp , list( stud.p ) , sum , na.rm=TRUE )
	wgt <- a1[,-1]
	# ISOP Model
	res.isop <- fit.isop( freq.correct , wgt , conv=conv , 
			maxit=maxit , progress=progress )	
    f1 <- freq.correct	
	f1 <- res.isop$fX
	# ADISOP Model (start with the fitted isop solution)
	res.adisop <- fit.adisop( freq.correct=f1 , wgt , conv=conv , maxit=maxit , 
			epsilon=epsilon , progress=progress  )
	wgt1 <- ( wgt / colSums( wgt ) ) / ncol(wgt)
	res.adisop$fit <- sqrt( sum( (  freq.correct - res.adisop$fX  )^2 * wgt1  ) )
	# logistic function
	res.logistic <- fit.logistic( freq.correct , wgt , 
						scores= scores / max( scores) , item.p=item.ps$x ,
						conv , maxit , progress=progress )
	# fit generalized logistic function
#	res.genlogistic <- fit.genlogistic( freq.correct , wgt , scores=scores / max( scores) , item.p=item.ps$x ,
#						conv , maxit , progress=progress )	
	# fit alpha1 parameter
	#*******
	# collect likelihood values
	ll <- data.frame( "model" = c("saturated" , "isop" , "adisop" , "logistic") , 
	   "ll" = c(res.isop$ll$ll.ind , res.isop$ll$ll.isop , 
				res.adisop$ll$ll.adisop , res.logistic$ll$ll.logistic ) )
	ll$llcase <- c( res.isop$ll$llcase.ind , res.isop$ll$llcase.isop , 
				      res.adisop$ll$llcase.adisop , 
					  res.logistic$ll$llcase.logistic  )
    #*****
    # model fit
	fit <- c( res.isop$fit , res.adisop$fit , res.logistic$fit )
	names(fit) <- c("isop","adisop","logistic")
	# probability tables
	prob.isop <- calc.prob.from.freq.dich( fittedfreq=res.isop$fX )	 
	prob.adisop <- calc.prob.from.freq.dich( fittedfreq=res.adisop$fX )	 
	prob.logistic <- calc.prob.from.freq.dich( fittedfreq=res.logistic$fX )	 
	prob.saturated <- calc.prob.from.freq.dich( fittedfreq= freq.correct )	 	

	#*******
	#* collect results
	res <- list( "freq.correct" = freq.correct , "wgt"=wgt ,
	          "prob.saturated" = prob.saturated , 
			  "prob.isop" = prob.isop , "prob.adisop" = prob.adisop ,
			  "prob.logistic" = prob.logistic,
			  "ll" = ll , "fit"=fit ,
			  "person"=person , "item"=item , 
			  "p.itemcat"=p.itemcat , "score.itemcat"=score.itemcat , 
			  "fit.isop" = res.isop , "fit.adisop" = res.adisop ,
			  "fit.logistic" = res.logistic , "s2" = Sys.time() , "dat" = dat ,
			  "model"="isop.dich")
	class(res) <- "isop"
	return(res)		
		}
##########################################################
# calculate probability table from fitted frequencies
calc.prob.from.freq.dich <- function( fittedfreq ){
		 SC <- nrow(fittedfreq)
		 I <- ncol(fittedfreq)
		 prob1 <- array( NA , dim=c(I,2,SC ) )
		 prob1[,2,] <- t(fittedfreq)
		 prob1[,1,] <- 1-t(fittedfreq)
		 return(prob1)	 
				}
#############################################################