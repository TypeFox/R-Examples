
###############################################
# nested multiple imputation
mice.nmi <- function( datlist , type="mice" , ... ){
    
	CALL <- match.call()
	ND <- length(datlist)
	imp <- list()
    
	cat("------ NESTED MULTIPLE IMPUTATION ------\n")
	for (dd in 1:ND){		
		data <- datlist[[dd]]
		cat("\n**************************************************\n")
		cat("***** Imputation Between Dataset" , dd , "******\n\n" )
		if ( type == "mice" ){
			imp[[dd]] <- mice::mice( data , ... )										
							}
		if ( type == "mice.1chain" ){
			imp[[dd]] <- mice.1chain( data , ... )										
							}							
			    }
	res <- list("imp" = imp , "type" = type )	
	if ( type=="mice"){
		NW <- imp[[1]]$m
					}
	if ( type=="mice.1chain"){
		NW <- imp[[1]]$Nimp
					}					
	res$imp[[1]]$call <- CALL				
	res$Nimp <- c ("between"=ND , "within" = NW )	
	res$nobs <- nrow(datlist[[1]])
	res$nvars <- ncol(datlist[[1]])
    class(res) <- "mids.nmi"			
	return(res)
		}
##################################################		
# summary method
summary.mids.nmi <- function( object , ... ){
	cat("Nested Multiple Imputation\n\n")
	Nimp <- object$Nimp
	cat( paste0("Number of between imputed datasets = " , Nimp["between"] , "\n") )
	cat( paste0("Number of within imputed datasets = " , Nimp["within"] , "\n") )	
    cat( paste0( object$nobs , " cases and " , object$nvars , " variables \n" ) )
	cat( paste0( "Number of iterations = " , object$imp[[1]]$iteration ) , "\n\n")
	imp0 <- object$imp[[1]]
	summary(imp0)
			}
####################################################			
# print method
print.mids.nmi <- function( x , ...){
	summary.mids.nmi( object=x , ...)
	}
#####################################################