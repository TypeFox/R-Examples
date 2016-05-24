##################################################
# attach all elements of an object in an environment   
.attach.environment.sirt <- function( res , envir ){
	CC <- base::length(res)
	for (cc in 1:CC){
		base::assign( base::names(res)[cc] , res[[cc]] , envir=envir )		
					}
			}
##################################################