
getInitial <- function(residual,phi,S){
     #.Call("getInitial", residual, phi, S, PACKAGE = "CNVLMM")
    temp <- .Call("getInitial", residual, phi, S)
#    cat('initial values: (sig2) ', temp[["sig2"]], ' (sig2g)', temp[["sig2g"]],'\n')
    unlist(temp)
}


doEM_REML <- function(curTheta,curK,y,X,yy,Xy,XX,phi,phiInv, trphiInv,S, q, threshold){
    res = .Call("doEM_REML", curTheta, curK, y, X, yy, Xy, XX, phi, phiInv, trphiInv, S, q, threshold)
    return(res)
}

doAI_REML<-function(curTheta,curK,y,X,yy,Xy,XX,phi,phiInv, trphiInv,S, q, itrmax,threshold){
    res = .Call("doAI_REML", curTheta, curK, y, X, yy, Xy, XX, phi, phiInv, trphiInv, S, q, itrmax, threshold)
    return(res)
}

doPolyGenic<-function(envirX,snp,pheno,phi,H0, thresEM,itrmax=20,thresAI,phiInv){

  S  <- length(pheno)
  if(H0) X  <- envirX else
  X  <- cbind(envirX,snp)
  y  <- matrix(pheno,ncol=1)
  yy <- t(y)%*%y
  Xy <- t(X)%*%y
  XX <- t(X)%*%X
  trphiInv<- sum(diag(phiInv))
  q  <- ncol(X)
  residual <- y - X%*%matrix(solve(XX)%*%Xy,ncol=1)

  inis     <- getInitial(residual,phi,S)

  if(inis[2]==0) {
      return(NA)
  } else {
  ## MME method to find the initial values
      updates  <- doEM_REML(inis[2],inis[1]/inis[2],y,X,yy,Xy,XX,phi,phiInv, trphiInv, S, q, thresEM)
  }
######

      curTheta <- updates[1]
      curK     <- updates[2]
      curTheta <- inis[2]
      curK     <- inis[1]/inis[2]

      cat('AI algorithm for REML.\n')
      results <- doAI_REML(curTheta,curK,y,X,yy,Xy,XX,phi,phiInv, trphiInv,S, q,itrmax,thresAI)

      return(results)
}





