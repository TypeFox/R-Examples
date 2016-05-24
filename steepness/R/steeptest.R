# Function to estimate p values for the steepness index based on the Dij or Pij measures #

steeptest <- function(X,rep,names=NULL,method=c("Dij","Pij"),order=TRUE){

# Is matrix X square? #

  if (nrow(X) != ncol(X))
    return("Error: Matrix X is not square and can not be analyzed")
  
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric");  
   
# Limit the number of replications #

  if ((rep < 1) | (rep > 1000000))
    return("Error: Number of replications must be between 1 and 1000000")

# Compute steepness measures for the original matrix #

  method <- match.arg(method)
  dyadc <- X + t(X);
  if (method == "Dij"){
    Dij <- X/dyadc-(((X/dyadc)-0.5)/(dyadc+1))
    Dij[is.nan(Dij)] <- 0.
    w1 <- rowSums(Dij);
    w2 <- Dij%*%w1;
    l1 <-colSums(Dij);
    l2 <- t(l1)%*%Dij;
  }
  if (method == "Pij"){
    Pij <- array(dim=c(nrow(X),ncol(X)),0.);
    Pij <- X/dyadc;
    Pij[is.nan(Pij)] <- 0.
    w1 <- rowSums(Pij);
    w2 <- Pij%*%w1;
    l1 <-colSums(Pij);
    l2 <- t(l1)%*%Pij;
  }
  DS <- w1 + w2 - l1 - t(l2);
  maxDS <- nrow(X)*(nrow(X)-1)/2;
  NormDS <- (DS + maxDS)/nrow(X);
  SortNormDS <- sort(NormDS,decreasing=TRUE,index.return=TRUE)
  rnk <- 1:nrow(X)
  Stp <- abs(lm(SortNormDS$x ~ rnk)$coefficients[2])
  names(Stp)<-NULL
  interc <- lm(SortNormDS$x ~ rnk)$coefficients[1]
  names(interc)<-NULL

# Carrying out the statistical test by means of a C program #

  vecX <- c(t(X))

  if (method == "Dij")
  out <- .C("steep",as.double(vecX),
          as.integer(nrow(X)),
          as.integer(rep),
          res1=double(rep),
          PACKAGE="steepness")
          
  if (method == "Pij")
  out <- .C("steep2",as.double(vecX),
          as.integer(nrow(X)),
          as.integer(rep),
          res1=double(rep),
          PACKAGE="steepness")

  Stpsim <- out$res1;

  if (is.null(names)) names <- paste('Ind.',1:nrow(X))
    else names <- names

  if (order == TRUE){
    names <- names[SortNormDS$ix]
    DS <- array(DS[SortNormDS$ix],dim=c(nrow(X),1))
    NormDS <- array(SortNormDS$x,dim=c(nrow(X),1))
  }

  if (method == "Dij")
    dimnames(Dij) <- list(c(names),c(names))
  if (method == "Pij")
    dimnames(Pij) <- list(c(names),c(names))
  
  dimnames(DS) <- list(c(names),"DS Values")
  dimnames(NormDS) <- list(c(names),"NormDS Values") 

  if (method == "Dij") matdom <- Dij
    else matdom <- Pij

  Stpresults <- list(call=match.call(),names = names,
  rep=rep,method=method,matdom=matdom,DS=DS,
  NormDS=NormDS,Stp=Stp,interc=interc,Stpsim=Stpsim)
  class(Stpresults) <- "steeptest"
  Stpresults
}
