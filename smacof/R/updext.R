# update function for transformation of external variables
updext <- function(x,w,external,extvars,constraint){
  m <- ncol(external)
  # reconstruct weigh matrix C
  if (constraint == "linear"){
    C <- solve(crossprod(external,w%*%external),crossprod(external,w%*%x))
  } else if (constraint == "diagonal") {
    C <- diag(colSums(external*(w%*%x))/colSums(external*(w%*%external)))
  }
  # For updating external[,s] we need a value larger than the largest eigenvalue
  # of V kronecker CC'.  
  svdC <- svd(C)
  v2 <- 2*max(diag(w))* svdC$d[1]^2                     # 2*max(diag(w)) is an upperbound of the largest eigenvalue of w
  z <- external - (1/v2)*(w %*% (external %*% C - x)) %*% t(C) # Compute the unconstrained majorization update
  external.old <- external
  iord.prim <- list()
  for (s in 1:ncol(external)){
    #     tt <- transform(z[,s], extvars[[s]], normq = 0)     # Compute update for external variable s
    #     external[,s] <- tt$res*(n/sum(tt$res^2))^.5         # Make the external variable of length n
    tt.plus <- transform(z[,s], extvars[[s]], normq = 0)     # Compute update for external variable s
    tt.min  <- transform(-z[,s], extvars[[s]], normq = 0)    # Compute update for external variable s
    if (sum((tt.plus$res - z[,s])^2) < sum(((tt.min$res + z[,s]))^2) ) {
      #external[, s] <- tt.plus$res*(n/sum(tt.plus$res^2))^.5
      iord.prim[[s]] <- tt.plus$iord.prim                    # Retain the ordening if primary approach to ties        
    } else {
      #external[, s] <- tt.min$res*(n/sum(tt.min$res^2))^.5
      iord.prim[[s]] <- tt.min$iord.prim                     # Retain the ordening if primary approach to ties        
    }      
  }
  return(list(external = external, iord.prim = iord.prim))
}
