

##########################################################
# fixed loadings in tam.mml.2pl
#  B.fixed 	
#  An optional matrix with four columns for fixing B matrix entries 
#   in 2PL estimation. 1st column: item index, 2nd column: category, 
#   3rd column: dimension, 4th column: fixed value. 
tamaanify.proc.loadings.B.fixed <- function(res){
	lavpartable <- res$lavpartable
	Q <- res$Q
	maxcat <- res$maxcat
	items <- res$items
	B.fixed <- NULL
	res$est.variance <- FALSE
    ind <- which( ( lavpartable$op == "=~" ) & 
				  ( lavpartable$free == 0 ) )
	if ( length(ind) > 0 ){			  
		lav1 <- lavpartable[ ind , ]
		N1 <- nrow(lav1)  
		for (vv in 1:N1){
		#	vv <- 1
			lav1.vv <- lav1[vv,]
			i2 <- which( colnames(res$resp) == paste(lav1.vv$rhs) )
			i3 <- items[ paste(items$item) == paste(lav1.vv$rhs) , "ncat"] - 1
			i4 <- which( colnames(Q) == paste(lav1.vv$lhs) )
			B1 <- cbind( i2 , seq(1,i3) + 1 , i4 , seq(1,i3)*lav1.vv$ustart )
			B.fixed <- rbind( B.fixed , B1 )
				}
		# res$est.variance <- TRUE
		facs <- colnames(Q)
		lav1 <- lavpartable[ ( paste0(lavpartable$lhs) %in% facs ) & 
					( paste0(lavpartable$rhs) %in% facs )	&
					( paste0(lavpartable$op) == "~~" )	, ]
		lav1 <- lav1[ lav1$free > 0 , ]
		if ( nrow(lav1) > 0 ){ 
				res$est.variance <- TRUE
						}		
			
#	  colnames(B.fixed) <- c("item_index" , "cat" , "dim" , "value") 					
#	  rownames(B.fixed) <- items[ B.fixed[,"item_index"] ]
	  
				}
	  res$B.fixed <- B.fixed	
	  
	  
	  return(res)
				}
#############################################################

