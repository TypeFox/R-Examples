#################################################
# create datlist
datlist_create <- function(datasets){
	CALL <- match.call()
	# if ( class(datasets) %in% c("mids","mids.1chain") ){
	if ( inherits( datasets , "mids")  |  inherits( datasets , "mids.1chain")
					){ 	
		datasets <- mids2datlist(datasets)
						}
	# if ( class(datasets) %in% "imputationList" ){
	if ( inherits(datasets , "imputationList" ) ) {
		datasets <- datasets$imputations								
						}
    class(datasets) <- "datlist"
	attr(datasets,"Nimp") <- length(datasets)
	attr(datasets,"call") <- CALL
	attr(datasets,"nobs") <- nrow(datasets[[1]])
	attr(datasets,"nvars") <- ncol(datasets[[1]])	
    return(datasets)
                    }
#**************** print method ***********************			
print.datlist <- function(x,...){
  cat("Object of class 'datlist'\nCall: ")
  print( attr(x,"call"))  
  cat("MI data with", attr(x,"Nimp") ,"datasets\n")
  v1 <- paste0( attr(x,"nobs") , " cases and " ,
	attr(x,"nvars") , " variables \n" )
  cat(v1)
}
########################################################
					

#######################################################
# create nested.datlist
nested.datlist_create <- function(datasets){
    CALL <- match.call()
	if ( class(datasets) %in% c("mids.nmi") ){
		datasets <- mids2datlist(datasets)
						}
	if ( class(datasets) %in% "NestedImputationList" ){
		datasets <- datasets$imputations								
						}
	v1 <- c("between" = length(datasets) , "within" = length(datasets[[1]]) )				
    class(datasets) <- "nested.datlist"
	attr(datasets,"Nimp") <- v1 			
	attr(datasets,"call") <- CALL
	attr(datasets,"nobs") <- nrow(datasets[[1]][[1]])
	attr(datasets,"nvars") <- ncol(datasets[[1]][[1]])			
    return(datasets)
                    }
#**************** print method ***********************					
print.nested.datlist <- function(x,...){ 
  cat("Object of class 'nested.datlist'\nCall: ")
  print( attr(x,"call"))  
  v1 <- paste0( "NMI data with ", attr(x,"Nimp")[1] ," between datasets and " ,
		  attr(x,"Nimp")[2] , " within datasets\n")
  cat(v1)    
  v1 <- paste0( attr(x,"nobs") , " cases and " ,
	attr(x,"nvars") , " variables \n" )
  cat(v1)  
}
#######################################################
					