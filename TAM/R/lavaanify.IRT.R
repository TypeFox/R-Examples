
########################################################################
# lavaanify extension
lavaanify.IRT <- function( lavmodel , items=NULL , data = NULL , include.residuals=TRUE ,
          doparse=TRUE ){
 a0 <- Sys.time()
 
    #*** check for input errors of data
	ind <- is.matrix( items ) | is.data.frame( items )
    if ( ind ){
	    items <- colnames(items)
				}

	#***
	if ( doparse ){
	   lavmodel <- doparse( lavmodel )
					}
 
    if ( is.null(items) ){
             items <- colnames(data)
					}	
	# grep for MEASERR
	lavmodel <- lavaanify.grep.MEASERR( lavmodel )										
	# grep for nonlinear terms
	res <- lavaanify.grep.nonlinear( lavmodel , items )
	lavmodel <- res$lavmodel
	nonlin_factors <- res$nonlin_factors
	nonlin_syntable <- res$nonlin_syntable
	items <- res$items
    res <- lavaanify.sirt.v1( lavmodel = lavmodel )	
	
	ind <- grep( "__[A-Z,a-z]" , lavmodel )
	
	if ( length(ind) > 0 ){
		res <- lavpartable.grep.underbrace( lavpartable=res$lavpartable , items )		
		
 		res <- remove.duplicated.variances.lavsyn(res , items)	
		res <- lavaanify.sirt.v1( lavmodel = res)
# cat("\n*** sirt.v1 second") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1					
		lavpar <- res$lavpartable
		lavsyn <- res$lavaan.syntax

		lavpar0 <- lavpar
		ind1 <- which( paste(lavpar$op) == "?="  )
		ind2 <- which( !( paste(lavpar$rhs) %in% c("g1","s1")  ) )
		ind <- intersect( ind1 , ind2 )
		if ( length(ind) > 0 ){
		  lavpar0 <- lavpar0[ - ind , ]
						}
		res$lavpartable <- lavpar0
		lavsyn1 <- unlist( strsplit( lavsyn , "\n") )
		cn <- items
		v2 <- paste0( cn , "?=1*" , cn )
		v2 <- c(v2,paste0( cn , "?=0*" , cn ))
		lavsyn1 <- lavsyn1[ ! (lavsyn1 %in% v2 ) ]
		lavsyn1 <- paste0( lavsyn1 , collapse="\n" )
		res$lavaan.syntax <- lavsyn1
 # cat("rest ind") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							
				}
		

	#****************
	# estimate residual variances
	if (include.residuals){
		lavpartable <- res$lavpartable
		ind <- which( ( paste(lavpartable$lhs) == paste(lavpartable$rhs) ) & 
				( paste(lavpartable$op) == "~~" ) & ( paste(lavpartable$lhs) %in% items ) )
		ind2 <- which( ( lavpartable$free == 0 ) & (lavpartable$ustart == 0 ) )
		ind <- intersect( ind , ind2 )
		if ( length(ind) > 0 ){
			LI <- length(ind)
			lavpartable[ind , "ustart"] <- NA
			lavpartable[ind , "free"] <- 1:LI  + 1000
							}
		lavsyn1 <- lavpartable2lavsyntax( lavpartable )								
		res$lavpartable <- lavpartable
		lavsyn1 <- paste0( lavsyn1 , collapse="\n" )
		res$lavaan.syntax <- lavsyn1

				}				
	# eliminate some entries of "?=" from parameter table
	# e.g. I1 ?= 1*I1
	lavpar0 <- res$lavpartable
	ind <- which( ( lavpar0$lhs == lavpar0$rhs ) & ( lavpar0$op == "?=" ) )
	if ( length(ind) > 0 ){	lavpar0 <- lavpar0[ -ind, ] }
	res$lavpartable <- lavpar0
	res$nonlin_factors <- nonlin_factors
	res$nonlin_syntable <- nonlin_syntable
# cat("all out") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1					
	return(res)
			}

			
##*****************************************************************			
###################################################################


