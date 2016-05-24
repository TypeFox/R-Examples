unitsimpls = function(Xtrain, Ytrain, Xtest = NULL, ncomp = NULL, weight = FALSE, beta) 
{
  pls.out <- simpls(Xtrain = Xtrain, Ytrain = Ytrain, Xtest = Xtest, ncomp, weight = weight, beta = beta)
  euclidian.norm <- function(xvec) {
      return(sqrt(sum(xvec * xvec)))
  }
  R.norm <- apply(pls.out$R, 2, euclidian.norm)
  if (length(R.norm) == 1) {
      M <- matrix(R.norm, 1, 1)
      Mi <- matrix(1/R.norm, 1, 1)
  }
  else {
      M <- diag(R.norm)
      Mi <- diag(1/R.norm)
  }
  Rnew <- pls.out$R %*% Mi
  Tnew <- pls.out$T %*% Mi
  Qnew <- pls.out$Q %*% M
  Pnew <- pls.out$P %*% M
  if (!is.null(Xtest)) 
      if (weight == TRUE){
          list(B = pls.out$B, Ypred = pls.out$Ypred, P = Pnew, 
               Q = Qnew, T = Tnew, R = Rnew, meanX = pls.out$meanX, weight = pls.out$weight)
      }
      else { 
          list(B = pls.out$B, Ypred = pls.out$Ypred, P = Pnew, 
               Q = Qnew, T = Tnew, R = Rnew, meanX = pls.out$meanX)
      }    
  else 
      if (weight == TRUE){
          list(B = pls.out$B, P = Pnew, Q = Qnew, T = Tnew, R = Rnew, 
               meanX = pls.out$meanX, meanY = pls.out$meanY, weight = weight)
      }
      else {
          list(B = pls.out$B, P = Pnew, Q = Qnew, T = Tnew, R = Rnew, 
               meanX = pls.out$meanX, meanY = pls.out$meanY)
      } 
}
