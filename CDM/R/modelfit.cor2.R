


###############################################################################
## extern "C" {
## modelfit_cor2__C( posterior_, data_, dataresp_, probs1_, probs0_, ip_, expiijj_) ;
## }
## modelfit_cor2_Cpp
modelfit_cor2_aux <- function( posterior_ , data_, dataresp_, probs1_ , 
			probs0_, ip_, expiijj_ ){
    # note that ip_ = ip -1 ;
	res <- .Call("modelfit_cor2_Cpp", 
					as.matrix(posterior_) , as.matrix(data_) , 
					as.matrix(dataresp_) , as.matrix(probs1_) , 
					as.matrix(probs0_) , as.matrix(ip_) , as.matrix(expiijj_) ,
					PACKAGE = "CDM")
    return(res)			
			}
###############################################################################


#############################################################################
modelfit.cor2 <-
function( data , posterior , probs ){
	K <- max( apply( data , 2 , max , na.rm=TRUE ) )
	if ( K>1 ){ stop("modelfit.cor only allows for dichotomous data\n") }
    data.resp <- 1 - is.na(data)
    data[ is.na(data) ] <- 9
    data1 <- data*data.resp
    I <- ncol(data)
    # calculate counts (ignore weights here!!) 
    n11 <- t(  ( data==1) * data.resp ) %*% ( ( data==1) * data.resp )
    n10 <- t(  ( data==1) * data.resp ) %*% ( ( data==0) * data.resp )
    n01 <- t(  ( data==0) * data.resp ) %*% ( ( data==1) * data.resp )
    n00 <- t(  ( data==0) * data.resp ) %*% ( ( data==0) * data.resp )
    
    #********************************
    # covariances 
    # library combinat is needed here
    ip <- itempairs <- t( utils::combn(I,2 ) )
    colnames(itempairs) <- c("item1" , "item2" )
    itempairs <- as.data.frame( itempairs )            
    itempairs$n11 <- n11[ ip ]
    itempairs$n10 <- n10[ ip ]
    itempairs$n01 <- n01[ ip ]
    itempairs$n00 <- n00[ ip ]
	itempairs$n <- rowSums( itempairs[ , c("n11","n10", "n01","n00") ] )
 #   itempairs$Exp00 <- itempairs$Exp01 <- itempairs$Exp10 <- itempairs$Exp11 <- NA
 #   itempairs$corExp <- itempairs$corObs <- NA
 
	#****
	# calculate expected score for every person and every item
	exp.ii.jj <- posterior %*% t( probs[,2,] )

	probs1 <- probs[ , 2 , ]
	probs0 <- probs[ , 1 , ]
    res <- modelfit_cor2_aux( posterior_ = posterior , data_ = data , 
			dataresp_ = data.resp , probs1_ =probs1 , probs0_ = probs0 , 
			ip_ = ip-1 , expiijj_ = exp.ii.jj  )
	r1 <- res$itempair_stat
	

	itempairs$Exp11 <- r1[,1]
	itempairs$Exp10 <- r1[,2]	
	itempairs$Exp01 <- r1[,3]
	itempairs$Exp00 <- r1[,4]	


	# observed correlation
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
#	itempairs$G2 <- NA	
	itempairs$RESIDCOV <- NA		
	itempairs$Q3 <- res$Q3			
	    
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
	
	itempairs <- itempairs[ itempairs$n > 0 , ]
	
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
			mean( abs( itempairs$corObs - itempairs$corExp ) , na.rm=TRUE) ,
			sqrt( mean( ( itempairs$corObs - itempairs$corExp )^2 , na.rm=TRUE ) ) ,			
			mean( itempairs$X2 ) , # mean( itempairs$G2) ,
			mean( 100*abs(itempairs$RESIDCOV ) , na.rm=TRUE ) ,
			mean( abs( itempairs$Q3 ) , na.rm=TRUE) ,
			mean( abs( itempairs$Q3 - mean(itempairs$Q3,na.rm=TRUE) ) , na.rm=TRUE )
						) )
	rownames(modelfit) <- c("MADcor" , "SRMSR" , "MX2" , # "MG2",
				"100*MADRESIDCOV" , "MADQ3" , "MADaQ3" )
    
#	modelfit <- modelfit[ ! ( rownames(modelfit) %in% 
#				c("MX2" , "MADQ3" , "MADaQ3" ) ) , , drop=FALSE ]
	modelfit <- modelfit[ ! ( rownames(modelfit) %in% 
				c("MX2"  ) ) , , drop=FALSE ]

	
	#*****
	# summary statistics
	modelfit.test <- data.frame("type" = c("max(X2)","abs(fcor)") , 
			"value" = c( max( itempairs$X2) , max( abs(itempairs$fcor) )  ) ,
			"p" = c( min( itempairs$X2_p.holm) , min( itempairs$fcor_p.holm)  ) 
				)
	#**** statistics for use in IRT.compareModels
	statlist <- data.frame( "maxX2"= modelfit.test[1,"value"] ,
					"p_maxX2"= modelfit.test[1,"p"] )	
	h1 <- modelfit$est	
	statlist <- cbind( statlist , t(h1 ) )
	names(statlist)[-c(1:2) ] <- rownames(modelfit)
	#****
	# print results
    res <- list( "modelfit.stat" = modelfit , "itempairs" = itempairs , 
		"modelfit.test" = modelfit.test  , "statlist" = statlist )
    return(res)
    }

#######################################################################	
	