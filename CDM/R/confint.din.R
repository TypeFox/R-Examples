
#############################################################
# confidence interval for objects of class din
confint.din <- function(object, parm , level=.95 , extended=FALSE, ind.item.skillprobs=TRUE, 
       ind.item= FALSE, diagcov = FALSE , h=.001 ,...){

	if ( ! missing( parm) ){ 	
		partable <- object$partable
		g1 <- parm %in% partable$parnames 
		if ( mean(g1) < 1 ){ 
			stop("Not all requested parameters in parameter table!\n")
							}
		p1 <- partable[ partable$free  ,]
		g1 <- parm %in% p1[,"parnames"] 
		if ( mean(g1) < 1 ){ 
			extended <- TRUE
							}
						}                 
   v1 <- vcov( object, extended=extended, infomat=FALSE ,ind.item.skillprobs=ind.item.skillprobs, 
                   ind.item= ind.item , diagcov = diagcov , h=h ,...)
   c1 <- attr( v1 , "coef" )                
   c1 <- c1[ ! duplicated(names(c1) ) ]
   se1 <- sqrt( diag( v1 ) )
   q1 <- stats::qnorm( 1 - (1-level)/2 )
   res <- data.frame( c1 - q1*se1 , c1+q1*se1)
   rownames(res) <- names(c1)
   colnames(res) <- c("2.5 %" , "97.5 %")
   if ( ! missing(parm) ){
        res <- res[ parm , ]
                        }  
   return(res)
        }    