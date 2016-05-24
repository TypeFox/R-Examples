
#######################################################
# list of nested of multiply imputed datasets
NestedImputationList <- function(datasets){
      CALL <- match.call()
	  Nimp <- c( length(datasets) , length(datasets[[1]] ) )	
	  names(Nimp) <- c("Between" , "Within")
	  rval<-list(imputations=datasets, "Nimp"=Nimp)
	  class(rval)<-"NestedImputationList"
	  rval$call <- CALL
  	  attr(rval,"nobs") <- nrow(datasets[[1]][[1]])
	  attr(rval,"nvars") <- ncol(datasets[[1]][[1]])		  
	  return(rval)
	}
#################################################################
#**************** print method ***********************					
print.NestedImputationList <- function(x,...){ 
  cat("Object of class 'NestedImputationList'\nCall: ")
  print( x$call )  
  v1 <- paste0( "NMI data with ", x$Nimp[1] ," between datasets and " ,
		  x$Nimp[2] , " within datasets\n")
  cat(v1)    
  v1 <- paste0( attr(x,"nobs") , " cases and " ,
	attr(x,"nvars") , " variables \n" )
  cat(v1)
	}
#####################################################