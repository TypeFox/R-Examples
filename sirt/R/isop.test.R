
############################################################
# testing the ISOP model
isop.test <- 
function( data , jackunits = 20 , weights=rep(1,nrow(data)) ){
	dat <- as.matrix( data )
	dat.resp <- 1 - is.na( dat )
	dat[ is.na(dat) ] <- 0
	weights <- rep( 1 , nrow(dat) )
	N <- nrow(dat) 
	if ( length(jackunits==1)){
		jackunits <- ( 1:N ) %% jackunits
							}
	JJ <- max( jackunits )+1    

	# apply Cpp routine
	res <- isop_tests_cpp( dat , dat.resp , weights , jackunits , JJ )

	#**********************
	# arrange output
	I <- ncol(dat)
	# Jackknife estimates
	itemstat <- as.data.frame( matrix( NA , nrow=I+1 , ncol=3 ) )
	colnames(itemstat) <- c("parm" , "est" , "se" )
	W1test <- res$W1test
	itemstat[1,1] <- "test"
	itemstat[1,2] <- W1test[1]
	itemstat[1,3] <- sqrt( sum( ( W1test[-1] - mean( W1test[-1] ) )^2 ) * (JJ-1) / JJ )
	W1i <- res$W1i
	itemstat[-1,1] <- colnames(dat)
	itemstat[-1,2] <- W1i[,1]
	itemstat[-1,3] <- sqrt( rowSums( ( W1i[,-1] - rowMeans( W1i[,-1] ) )^2 ) * (JJ-1) / JJ )
	# include N's and p values 
	itemstat <- data.frame("parm" = itemstat$parm , "N" = c(N , colSums( dat.resp * weights ) ) , 
			"M" = c(NA , colSums( dat.resp * weights * dat ) / colSums( dat.resp * weights ) ) ,
			itemstat[,-1] )
	itemstat$t <- itemstat$est / itemstat$se        
	rownames(itemstat) <- NULL
	res <- list("itemstat" = itemstat , "Es" = res$Esi[,1], "Ed" = res$Edi[,1])
	res$JJ <- JJ
	class(res) <- "isop.test"
	return(res)
}
####################################################
# summary for ISOP test
summary.isop.test <- function( object , ... ){
    obji <- object$itemstat
	VV <- ncol(obji)
	cat("*** Test for the W1 Axiom in the ISOP Model **** \n\n")
	for (vv in 2:VV){
		obji[,vv] <- round( obji[,vv],3) 
				}
	print(obji)
	cat(paste0("\n-- Statistical inference is based on ", object$JJ ,
			" jackknife units.\n"))
	}
##########################################################
# call to Rcpp function
# res <- isop_tests_cpp( dat , dat.resp , weights , jackunits , JJ )
# SEXP isop_tests_C( SEXP dat_, SEXP dat_resp_, SEXP weights_, SEXP jackunits_, SEXP JJ_) ;
isop_tests_cpp <- function ( dat , dat.resp , weights , jackunits , JJ ){ 
		.Call("isop_tests_C", 
				 dat_ = dat ,  dat_resp_ = dat.resp ,  weights_ = weights ,  
				 jackunits_ = jackunits ,  JJ_=JJ ,
				PACKAGE = "sirt")
					}	
#############################################################