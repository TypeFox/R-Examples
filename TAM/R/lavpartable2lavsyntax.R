

#######################################################################
# create a lavaan syntax from lavaan parameter table
# This syntax only works for single groups.
lavpartable2lavsyntax <- function( lavpartable ){
	lavpartable <- lavpartable[ lavpartable$user != 2 , ]
    LL <- nrow(lavpartable)    
    syn0 <- paste0( lavpartable$lhs , lavpartable$op )	
	lavpartable$ustart2 <- lavpartable$ustart
	lavpartable$ustart2[ is.na( lavpartable$ustart) ] <- ""
	lavpartable$ind0 <- 0
	# op == "==" and ustart == NA
	ind <- which( ( lavpartable$op == "==" ) & ( is.na( lavpartable$ustart ) ) )	
	if ( length(ind) > 0 ){
		lavpartable[ ind , "ustart2" ] <- ""
		lavpartable[ ind , "ind0" ] <- 1
							}
			
	lavpartable$ustart2 <- 
			ifelse( lavpartable$ustart2 != "" , paste0(lavpartable$ustart2 , "*"  )	, lavpartable$ustart2 )						
														
    lavpartable$prefix <- ""		
    lavpartable$prefix <- ifelse( paste(lavpartable$label ) != "" , 
                paste0(lavpartable$label , "*"  ) , lavpartable$prefix ) 
    lavpartable$prefix <- ifelse( ( lavpartable$free == 0 ) & ( paste(lavpartable$label ) == "" ) , 
                 lavpartable$ustart2 , lavpartable$prefix )              
    syn0 <- paste0( syn0 , lavpartable$prefix , lavpartable$rhs  )
    syn0 <- paste0( syn0 , collapse="\n")
    return(syn0)
            }			
#######################################################################
	