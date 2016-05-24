
#########################################################################
IRT.se.din <- function( object , extended = FALSE , parm=NULL ,
            level = .95 , infomat=FALSE , ind.item.skillprobs = TRUE , 
			ind.item= FALSE , diagcov = FALSE , h=.001 , ... ){
			
	if ( ! missing( parm) ){ 	
		partable <- object$partable
		g1 <- parm %in% partable$parnames 
		if ( mean(g1) < 1 ){ 
			stop("Not all requested parameters in parameter table!\n")
							}
						}			
	# variance covariance matrix
    v1 <- vcov( object, extended=extended, infomat=FALSE ,ind.item.skillprobs=ind.item.skillprobs, 
                   ind.item= ind.item , diagcov = diagcov , h=h ,...)			
	# confidence intervals
	res <- data.frame( "parameter" = names(attr( v1 , "coef")) )
	res$est <- attr(v1,"coef")		
	res$se <- sqrt( diag(v1)  )
	
	level1 <-  round( 100*( 1  - level ) / 2 , 1 )
	quant1 <- abs( stats::qnorm( ( 1-level) / 2) )
	res[ , paste0( level1 , " %" ) ] <- res$est - quant1*res$se
	res[ , paste0( 100-level1 , " %" ) ] <- res$est + quant1*res$se	
	
	partable <- object$partable	
	ind1 <- match( res$parameter , partable$parnames )	
	vars1 <- c("partype" , "parindex" )
	res <- cbind( partable[ ind1, vars1 ] , res )
	vars2 <- c( "item" , "item.name" , "skillclass" , "fixed" ,"free" , 
				"rule" , "totindex"  )
	res <- cbind( res , partable[ ind1, vars2 ])	
    return(res) 
			}
			
			