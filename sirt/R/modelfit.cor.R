
##############################################################
# summary
summary.modelfit.sirt <- function( object , ... ){	
	cat("Test of Global Model Fit\n")
	obji <- object$modelfit.test
	for (vv in seq(2,ncol(obji))){	obji[,vv] <- round( obji[,vv] , 5 ) }
	print(obji)
	cat("\nFit Statistics\n")
	obji <- object$modelfit.stat
	for (vv in seq(1,ncol(obji))){	obji[,vv] <- round( obji[,vv] , 5 ) }
	print(obji)		
		}
#################################################################	



##########################################
# Modelfit in sirt
modelfit.sirt <- function( object ){
	#****
	# object of class tam.mml, tam.mml.2pl or tam.fa
	# Note that only dichotomous responses are allowed
	if (class(object) %in% c("tam.mml","tam.mml.2pl")){
		mod <- object
		# fmod1b <- modelfit.cor2(data=dat, 
		#	posterior=mod1b$post, probs=mod1b$rprobs)
		posterior <- mod$hwt
		probs <- mod$rprobs
		dat <- mod$resp
		dat[ mod$resp.ind == 0 ] <- NA
					}
	#*****
	# rasch.mml
	if (class(object)=="rasch.mml"){
		mod <- object
		posterior <- mod$f.qk.yi
		prob1 <- mod$pjk
			probs <- array( NA , dim=c( ncol(prob1) , 2 , nrow(prob1)) )
			probs[ , 2 , ] <- t(prob1)
			probs[ , 1 , ] <- 1 - t(prob1)
		dat <- mod$dat
					}
	#*****
	# rasch.mirtlc
	if (class(object)=="rasch.mirtlc"){
		mod <- object$estep.res
		posterior <- mod$f.qk.yi
		prob1 <- mod$pjk
			probs <- array( NA , dim=c( ncol(prob1) , 2 , nrow(prob1)) )
			probs[ , 2 , ] <- t(prob1)
			probs[ , 1 , ] <- 1 - t(prob1)
		dat <- object$dat
					}	
	#******
	# rasch.pml
	if ( class(object) !="rasch.pml"){ pmlobject <- NULL } else {
		data <- NULL ; posterior <- NULL ; probs <- NULL ; pmlobject <- object }	
	#*******
	# smirt	
	if (class(object) == "smirt"){
		# note that for polytomous response data some adaptations are
		# necessary: see modelfit in the CDM package
		mod <- object
		probs <- mod$probs
		posterior <- mod$f.qk.yi
		dat <- mod$dat
					}
	#*******
	# smirt	
	if (class(object) == "gom"){
		mod <- object
		probs <- mod$probs
		posterior <- mod$f.qk.yi
		dat <- mod$dat
					}					
	#*******
	# rm.facets
#	if (class(object) %in% c("rm.facets") ){
#		mod <- object
#		probs <- mod$probs
#		posterior <- mod$f.qk.yi
#		dat <- mod$procdata$dat2.NA
#					}						
					
	#*******
	# mirt	
	if (class(object) == "ConfirmatoryClass" | class(object)=="ExploratoryClass" ){
		mod <- object
		mod <- mirt.wrapper.posterior(mod)		
		probs <- mod$probs
		posterior <- mod$f.qk.yi
		dat <- mod$dat
					}					
	#******
	# R2noharm, noharm.sirt
    if ( class(object) %in% c("R2noharm","noharm.sirt") ){
		# exclusion criteria for noharm.sirt
		if ( class(object) == "noharm.sirt") {
			if ( object$estpars$estPsi > 0 ){
				stop("Model fit cannot be calculated because of correlated residuals")
										}
			if ( ! ( object$wgtm.default ) ){
				stop("Model fit cannot be calculated because not all item pairs are used for estimation")
										}										
								}
		  # evaluation of posterior
		  mod <- R2noharm.EAP(noharmobj=object, theta.k = seq(-6, 6, len = 15 ) ,
				print.output=FALSE )
		  probs <- aperm( mod$probs , c(1,3,2) )
		  posterior <- mod$posterior
		  dat <- object$dat		  		  
									}
	# calculate modelfit.cor
    if ( class(object) == "rasch.pml" ){
		res <- modelfit.cor.sirt.pml( data = dat , posterior =posterior , probs = probs ,
				           pmlobject=pmlobject)
								} else {
		res <- CDM::modelfit.cor2( data = dat , posterior =posterior , probs = probs )								
							}
    class(res) <- "modelfit.sirt"							
	return(res)
	}
################################################################################

#############################################################################
modelfit.cor.sirt.pml <-
function( data=NULL , posterior=NULL , probs=NULL , pmlobject=NULL ){

		data <- pmlobject$data
		itempairs <- as.data.frame( pmlobject$itempairs )            
		ip0 <- itempairs
		itempairs$n11 <- itempairs$f11
		itempairs$n10 <- itempairs$f10
		itempairs$n01 <- itempairs$f01
		itempairs$n00 <- itempairs$f00
		itempairs$n <- rowSums( itempairs[ , c("n11","n10", "n01","n00") ] )		
		itempairs$Exp00 <- ip0$p00 * itempairs$n
		itempairs$Exp10 <- ip0$p10 * itempairs$n
		itempairs$Exp01 <- ip0$p01 * itempairs$n
		itempairs$Exp11 <- ip0$p11 * itempairs$n

		n <- itempairs$n
		m1 <- ( itempairs$n10 + itempairs$n11 ) / n
		m2 <- ( itempairs$n01 + itempairs$n11 ) / n
		t1 <- itempairs$n11 / n  - m1 * m2
		itempairs$corObs <- t1 / sqrt( m1 * ( 1 - m1 ) * m2 * ( 1-m2 ) )
	 
		# observed correlation
		m1 <- ( itempairs$Exp10 + itempairs$Exp11 ) / n
		m2 <- ( itempairs$Exp01 + itempairs$Exp11 ) / n
		t1 <- itempairs$Exp11 / n  - m1 * m2
		itempairs$corExp <- t1 / sqrt( m1 * ( 1 - m1 ) * m2 * ( 1-m2 ) ) 

	#    m1 <- matrix( c(1,1,1,0,0,1,0,0) , 4 , 2 , byrow=T )
		
		# define further quantities
		itempairs$X2 <- NA
		itempairs$RESIDCOV <- NA			
			
		##############################
		itempairs$X2 <- ( itempairs$n00 - itempairs$Exp00 )^2 / itempairs$Exp00 +
						( itempairs$n10 - itempairs$Exp10 )^2 / itempairs$Exp10	+
						( itempairs$n01 - itempairs$Exp01 )^2 / itempairs$Exp01	+
						( itempairs$n11 - itempairs$Exp11 )^2 / itempairs$Exp11						
		
		itempairs$RESIDCOV <- ( itempairs$n11 * itempairs$n00 - itempairs$n10 * itempairs$n01 ) / itempairs$n^2 -
					( itempairs$Exp11 * itempairs$Exp00 - itempairs$Exp10 * itempairs$Exp01 ) / itempairs$n^2	
		##############################
	
		# labels
		itempairs$item1 <- colnames(data)[ itempairs$item1 ]
		itempairs$item2 <- colnames(data)[ itempairs$item2 ]

		# fisherz from psych package
		# residual of correlation
		itempairs$fcor <- psych::fisherz( itempairs$corObs ) - psych::fisherz( itempairs$corExp )
		
		#----
		# p values and p value adjustments adjustments
		
		# X2 statistic
		itempairs$X2_df <- 1
		itempairs$X2_p <- 1 - stats::pchisq(itempairs$X2 , df=1 )
		itempairs$X2_p.holm <- stats::p.adjust( itempairs$X2_p , method="holm")
		itempairs$X2_sig.holm <- 1 * ( itempairs$X2_p.holm < .05 )	
		itempairs$X2_p.fdr <- stats::p.adjust( itempairs$X2_p , method="fdr")
		# fcor statistic
		itempairs$fcor_se <- ( itempairs$n - 3 )^(-1/2)
		itempairs$fcor_z <- itempairs$fcor / itempairs$fcor_se
		itempairs$fcor_p <- 1 - stats::pnorm( abs(itempairs$fcor_z ) )
		itempairs$fcor_p.holm <- stats::p.adjust( itempairs$fcor_p , method="holm")
		itempairs$fcor_p.fdr <- stats::p.adjust( itempairs$fcor_p , method="fdr")

		#**********************
		# model fit
		modelfit <- data.frame( "est" = c( 
				mean( abs( itempairs$corObs - itempairs$corExp ) ) ,
				sqrt( mean( ( itempairs$corObs - itempairs$corExp )^2 ) ) ,			
				mean( itempairs$X2 ) , # mean( itempairs$G2) ,
				mean( 100*abs(itempairs$RESIDCOV ) ) 
						)
							) 
		rownames(modelfit) <- c("MADcor" , "SRMSR" , "MX2" , # "MG2",
					"100*MADRESIDCOV" )
		
	#*****
	# summary statistics
	modelfit.test <- data.frame("type" = c("max(X2)","abs(fcor)") , 
			"value" = c( max( itempairs$X2) , max( abs(itempairs$fcor) )  ) ,
			"p" = c( min( itempairs$X2_p.holm) , min( itempairs$fcor_p.holm)  ) 
				)					
	#****
	# print results
    res <- list( "modelfit.stat" = modelfit , "itempairs" = itempairs , 
		"modelfit.test" = modelfit.test  )		
    return(res)
    }
#######################################################################

#######################################################################	
# auxiliary function: weighted correlation	
.corr.wt <- function( x, y, w = rep(1,length(x))) {
#  stopifnot(length(x) == dim(y)[2] )
  w <- w / sum(w)
  # Center x and y, using the weighted means
  x <- x - sum(x * w)
  ty <- y - sum( y * w)
  # Compute the variance
  vx <- sum(w * x * x)
  vy <- sum(w * ty * ty)
  # Compute the covariance
  vxy <- sum(ty * x * w)
  # Compute the correlation
  vxy / sqrt(vx * vy)
}
##################################################################################