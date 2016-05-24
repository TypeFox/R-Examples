### AVERAGE RELATIVE BIAS
## B = matrice parametri di bootstrap
## B0 true parameters
#rm(list=ls())
#load("~/lavori/Rdevel/testfile/arbtest.rda")

arb <- function(Bpar,B0) {
  Bpar <- na.omit(Bpar)
  B0den <- ifelse(B0==0,1,B0)
  
  if (is.null(ncol(Bpar))) {
    Bpar <- matrix(Bpar,ncol=1)
  } 

  P <- NULL
  for (j in 1:ncol(Bpar)) {
    P <- cbind(P,(Bpar[,j]-B0[j])/B0den[j])
  }
  P <- apply(P,1,mean,na.rm=TRUE)
  arb <- sum(P)*100/nrow(Bpar)    

  return(arb)
}

amse <- function(Bpar,B0) {
  Bpar <- na.omit(Bpar)
  B0den <- ifelse(B0==0,1,B0)
  
  if (is.null(ncol(Bpar))) {
    Bpar <- matrix(Bpar,ncol=1)
  } 
    
  P <- NULL
  for (j in 1:ncol(Bpar)) {
    P <- cbind(P,((Bpar[,j]-B0[j])/B0den[j])^2)
  }
  P <- sqrt(apply(P,1,mean,na.rm=TRUE))
  amse <- sum(P)/nrow(Bpar)    
  
  return(amse)
}

### modifica del 3/6/2014 aggiunta gestione dati di tipo non matrice.
