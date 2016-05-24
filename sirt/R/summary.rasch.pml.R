

#*******************************************************
# Summary for rasch.copula object                         *
##NS S3method(summary,rasch.pml)
summary.rasch.pml <- function( object , ...){
    cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    cat("  Function" , object$fct , " \n")
    cat("-----------------------------------------------------------------\n")
    cat("Pairwise Marginal Likelihood Estimation \n")
	cat(paste("Link function:" , object$link , "\n"))
    cat("-----------------------------------------------------------------\n")
    cat( "Pseudolikelihood objective function = " , round( object$ic$deviance , 2 ) , "\n" )
    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "Number of used item pairs in PML estimation = " , nrow(object$itempairs) , "\n" )  
#    cat( "AIC = " , round( object$ic$AIC , 2 ) , " ; penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , "\n" )    
    cat( "PLIC = " , round( object$ic$PLIC , 2 ) , " ; penalty =" , 
				round( object$ic$PLIC - object$ic$deviance ,2 ) , 
				" | PLIC = -2*PL + log(n)*p \n" )
    cat( "      (Pseudolikelihood information criterion)\n")				
#    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," ; penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) , "\n\n" )         
#    cat( "Trait Distribution (" , length(object$theta.k) , " Knots )\n" , 
#              "Mean=" , 0 , " SD=" , 1 , "\n") 

    cat("-----------------------------------------------------------------\n")
	# print item summary
	cat("Item Parameter Summary\n")
	cat( " Estimated" , length(object$bG)+length(object$aG) , "Item Parameters\n\n")
	.pr( object$item , digits=3 )		# print item statistics
	#....................................................................
	# print Trait parameter summary
    cat("-----------------------------------------------------------------\n")	
	cat("Trait SD (Probit Link): ")
    cat( round(object$sigma , 3 ) , "\n")
	cat("Trait SD (Logit Link) : ")
    cat( round( object$item$sigma[1] * 1.701 , 3 ) , "\n")		
	if( object$fct	== "rasch.pml3" ){
		cat("\nGreen-Yang Reliability \n")
		cat( "omega =" , round( object$omega.rel , 3 ) , "\n")
				}    
	
	cat("-----------------------------------------------------------------\n")
	
	if ( object$est.corrs){
    cat("-----------------------------------------------------------------\n")
	cat("Residual Correlation Parameter Summary\n")
		cat( " Estimated" , length(object$epsG) , "Residual Correlation Parameters\n\n")
	   .pr( object$eps.corrM , digits=3 )		# print item statistics
			}		
		}
#*******************************************************
