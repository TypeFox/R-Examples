regH <- function(H, type = 1){

  epsilon <- sqrt(.Machine$double.eps)

  if(type == 1){

      op <- options(warn = -1)
      R <- chol(H, pivot = TRUE)
      options(op)
      p <- dim(H)[2]
      ipiv <- piv <- attr(R, "pivot")
      ipiv[piv] <- 1:p
      rank <- attr(R, "rank")
      ind <- 1:rank
      if (rank < p) R[(rank + 1):p, ] <- 0
      R <- R[ipiv, ipiv]
      H <- crossprod(R)
  
} 


  if(type == 2){
  
  ds <- 0
  
  eH <- eigen(H, symmetric = TRUE)
  if(min(eH$values) < epsilon){ eH$values <- abs(eH$values); ds <- 1 }
  if(min(eH$values) < epsilon){ eH$values[which(eH$values < epsilon)] <- 0.0000001; ds <- 1 }
  
  if(ds == 1) H <- eH$vectors%*%tcrossprod(diag(1/eH$values),eH$vectors)   
  
  } 
  
  
  
  
  if(type == 3){
  
  D <- diag(abs(diag(H))^-.5)
  H1 <- D%*%H%*%D
  H1inv <- D%*%solve(H1)%*%D
  
  # does not solve much on a well specified model
  # we need to do this using the inverse
  
  
  #ds <- 0
  #
  #eH <- eigen(H, symmetric = TRUE)
  #if(min(eH$values) < epsilon){ eH$values <- abs(eH$values); ds <- 1 }
  #if(min(eH$values) < epsilon){ eH$values[which(eH$values < epsilon)] <- 0.0000001; ds <- 1 }
  #
  #if(ds == 1) H <- eH$vectors%*%tcrossprod(diag(1/eH$values),eH$vectors)   
  #
  }   
  
  
  H
     
}




     























