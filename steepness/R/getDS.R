# Function to obtain the David's scores (DS) based on dyadic dominance indices #

getDS <- function(X,names=NULL,method=c("Dij","Pij")){
  if (nrow(X) != ncol(X)) 
    return("Error: Sociomatrix must be square")
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric")
method <- match.arg(method)
dyadc <- X + t(X)
if (method == "Dij"){
  Dij <- array(dim=c(nrow(X),ncol(X)),0.)
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
DS <- w1 + w2 - l1 - t(l2);
if (is.null(names)) names <- paste('Ind.',1:nrow(X))
DS <- array(DS,dim=c(nrow(X),1),dimnames=c(list(names),"David's Scores"))
return(DS)
}

