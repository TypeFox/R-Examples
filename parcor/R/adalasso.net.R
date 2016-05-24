

`adalasso.net` <-
function(X, k=10,use.Gram=FALSE,both=TRUE,verbose=FALSE,intercept=TRUE)
{
  p <- ncol(X)
  
  X <- scale(X)     # data needs to be centered and standardized
  colnames(X) <- 1:p    # each column gets a name
                      
  B.lasso<-B.adalasso <- matrix(0, nrow=p, ncol=p)
  
  colnames(B.lasso) <- colnames(B.adalasso)<-1:p
  pcor.adalasso<-NULL
  if (verbose==TRUE){
  cat(paste("Performing local (adaptive) lasso regressions\n"))
  cat(paste("Vertex no "))
}
  for (i in 1:p) ## visit all nodes
  {
    if (verbose==TRUE){
    if ((i/10)==floor(i/10)) {
      cat(paste(i,"..."))}
  }
    noti <- (1:p)[-i]
    yi <- X[ ,i]       ## response
    Xi <- X[ ,noti]    ## predicted by all other nodes with i missing
    
    ## perform adaptive lasso regression & extract regression coefficients  
    dummy <- adalasso(Xi, yi, k=k,use.Gram=use.Gram,both=both,intercept=intercept)
    coefi.lasso<-dummy$coefficients.lasso
     B.lasso[i,-i] <- coefi.lasso
    if (both==TRUE){
    coefi.adalasso<-dummy$coefficients.adalasso
      B.adalasso[i,-i] <- coefi.adalasso 
    }   
    }
    pcor.lasso<-Beta2parcor(B.lasso,verbose=verbose)
    if (both==TRUE){
    pcor.adalasso <- Beta2parcor(B.adalasso,verbose=verbose)
  }
  cat(paste("\n"))
  

  return(list(pcor.lasso=pcor.lasso,pcor.adalasso=pcor.adalasso))  
}
