			
#################################################################
# grep for "__" operator, meaning I01__I10
lavpartable.grep.underbrace <- function( lavpartable , items ){
    lav2 <- lavpartable
	lav2 <- lav2[ lav2$user != 2 , ]
	LL <- nrow(lav2)
	
	lav2[ lav2$op == "==" , "free" ] <- -9
	
	syn <- NULL

    for (ll in 1:LL){
        # ll <- 13
        lav2.ll <- lav2[ll,]
		##*** ARb 2014-09-21: bug fix for "~1" operator
		if (lav2.ll$op == "~1" ){
			lav2.ll$rhs <- "1" 
			lav2.ll$op <- "~"
						}			
        ind.ll <- grep( "__" , c( lav2.ll$lhs  , lav2.ll$rhs  ) )
		
		#****** lhs
        v11 <- v1 <- lav2.ll$lhs            
        v10 <- strsplit( v1 , "__" )[[1]] 
        if ( length(v10) > 1 ){
              v11 <- items[ seq( which( items == v10[1] ) , 
									which( items == v10[2] )  ) ]            
                               }   
        syn0 <- paste0( v11 , " " , lav2.ll$op , " ")
        g1 <- ifelse( lav2.ll$label != "" ,  paste0( lav2.ll$label , "*" ) , "" )
        g1 <- paste0( g1 , "" ,  ifelse( lav2.ll$free == 0 ,  paste0( lav2.ll$ustart , "*" ) , "" ) )
        syn0 <- paste0( syn0 , g1 , "" )

        #******* rhs		
        v11 <- v1 <- lav2.ll$rhs            
        v10 <- strsplit( v1 , "__" )[[1]]
         if ( length(v10) > 1 ){        
                v11 <- items[ seq( which( items == v10[1] ) , which( items == v10[2] )  ) ]            
                                    }
        syn0 <- paste0( syn0 , paste0( v11 , "\n")    )
        syn0 <- paste0( syn0 , collapse = "" )		
        syn <- paste0( syn , syn0 , collapse = "")   
	
                }   

    return(syn)
        }
########################################################################
