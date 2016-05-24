
#####################################################################
MIcombine.NestedImputationResultList <- function(results, ...){

  Nimp <- c( length(results) , length( results[[1]]) )
  p1 <- coef(results[[1]][[1]])
  NV <- length( p1 )
 
  NB <- Nimp[1]
  NW <- Nimp[2]
  
  # collect estimated variance matrices
  qhat <- array( NA , dim=c( NB , NW , NV ) )
  dimnames(qhat)[[3]] <- names(p1) 
  dimnames(qhat)[[1]] <- paste0("Between_Imp" , 1:NB ) 
  dimnames(qhat)[[2]] <- paste0("Within_Imp" , 1:NW )  
  
  u <- array( NA , dim=c( NB , NW , NV , NV) )
  dimnames(u)[[4]] <- dimnames(u)[[3]] <- names(p1) 
  dimnames(u)[[1]] <- paste0("Between_Imp" , 1:NB ) 
  dimnames(u)[[2]] <- paste0("Within_Imp" , 1:NW )  
  for (bb in 1:NB){
    for (ww in 1:NW){
	     u[bb,ww,,] <- vcov( results[[bb]][[ww]] )
		 qhat[bb,ww,] <- coef( results[[bb]][[ww]] )
					}
				}
  rval <- pool.nmi.scalar.helper( qhat , u , NV , NB , NW )
  rval$Nimp <- Nimp
  # rval$vcov <- rval$Tm
  # rval$coef <- rval$qbar
  class(rval) <- "mipo.nmi"  
  return(rval) 
}
#################################################################