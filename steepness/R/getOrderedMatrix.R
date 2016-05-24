getOrderedMatrix <- function(X,names=NULL,method=c("Dij","Pij")){
  if (nrow(X) != ncol(X)) 
    return("Error: Sociomatrix must be square")
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric")
  method <- match.arg(method)
  dyadc <- X + t(X)
  if (method == "Dij"){
    Dij <- X/dyadc-(((X/dyadc)-0.5)/(dyadc+1))
    Dij[is.nan(Dij)] <- 0.
    w1 <- rowSums(Dij)
    w2 <- Dij%*%w1
    l1 <-colSums(Dij)
    l2 <- t(l1)%*%Dij
  }
  if (method == "Pij"){
    Pij <- array(dim=c(nrow(X),ncol(X)),0.)
    Pij <- X/dyadc
    Pij[is.nan(Pij)] <- 0.
    w1 <- rowSums(Pij)
    w2 <- Pij%*%w1
    l1 <-colSums(Pij)
    l2 <- t(l1)%*%Pij
  }  
  DS <- w1 + w2 - l1 - t(l2)
  maxDS <- nrow(X)*(nrow(X)-1)/2
  NormDS <- (DS + maxDS)/nrow(X)
  SortNormDS <- sort(NormDS,decreasing=TRUE,index.return=TRUE)
  matord <- X[SortNormDS$ix,SortNormDS$ix]
  namesord <- names[SortNormDS$ix]
  rownames(matord) <- namesord
  colnames(matord) <- namesord
  
  list(ordered.matrix=matord,ordered.names=namesord,order.seq=SortNormDS$ix)
  
}
  
  