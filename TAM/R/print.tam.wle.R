
#############################################
# print method for tam.wle objects
print.tam.wle <- function(x, digits = 3 , ...){
	ndim <- attr(x,"ndim")
	nobs <- attr(x,"nobs")
    CALL <- attr(x,"call")
	WLE.rel <- attr(x,"WLE.rel")
	#*** print general informations
	print_tam_wle_general( ndim , nobs, CALL )		
    #*** print reliabilities	
    print_tam_wle_WLErel( WLE.rel, ndim , digits )	
}
##############################################


#**********************************************************************
print_tam_wle_general <- function( ndim , nobs, CALL ){
	#*** print output with general informations
    cat("Object of class 'tam.wle'\nCall: ")
    print( CALL ) 
    # cat("\n")	
	if (ndim==1){ D0 <- ""} else { D0 <- "s" }
    v1 <- paste0("  WLEs for ", nobs , " observations and " ,
		    ndim , " dimension" , D0 , "\n")
	cat(v1)
		}
#**********************************************************************		
print_tam_wle_WLErel <- function( WLE.rel, ndim , digits ){		
	if (ndim==1){
		v1 <- paste0("  WLE Reliability = " , round(WLE.rel,digits=digits) , "\n")
		cat(v1)
		}
	if (ndim>1){
	   for (dd in 1:ndim){
		v1 <- paste0("  WLE Reliability (Dimension " , dd , 
					") = " , round(WLE.rel[dd],digits=digits) , "\n")
		cat(v1)	   
				}
		}		
	}
#**********************************************************************	