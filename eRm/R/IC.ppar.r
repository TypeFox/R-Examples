IC.ppar <- function(object)
{
#computes loglik, AIC, BIC, and cAIC  for JML, MML, CML
#object of class ppar


  #---------- full likelihood ----------
  X <- object$X
  if (length(object$pers.ex) > 0) X01 <- object$X01[-object$pers.ex,] else X01 <- object$X01
  mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
  mt_ind <- rep(1:length(mt_vek),mt_vek)
  mt_seq <- sequence(mt_vek)
  gmemb <- object$gmemb

  pmat <- pmat(object) 
  pmat.l0 <- tapply(1:length(mt_ind),mt_ind, function(ind) { #expand pmat for 0-th category
                             vec0 <- 1-rowSums(as.matrix(pmat[,ind]))     #prob for 0th category
                             cbind(vec0,pmat[,ind])
                            })
  pmat0 <- matrix(unlist(pmat.l0),nrow=length(gmemb))        #X01 matrix 0th category included
  X01.l0 <- tapply(1:length(mt_ind), mt_ind, function(ind) { #expand X01 for 0-th category
                            vec0 <- 1-rowSums(as.matrix(X01[,ind]))     #prob for 0th category
                            cbind(vec0,X01[,ind])
                            })
  X010 <- matrix(unlist(X01.l0),nrow=length(gmemb))          #X01 matrix 0th category included
  loglik.full <- sum(log(na.exclude(pmat0[X010 == 1])))      #vector of "observed" solving probabilities

  N.ex <- dim(object$X.ex)[1]                                #number of persons (excluded)
  npar.full <- (dim(object$W)[2])+sum(object$npar)           #number of item + person parameters
  AIC.full <- -2*loglik.full + 2*npar.full
  BIC.full <- -2*loglik.full + log(N.ex)*npar.full
  cAIC.full <- -2*loglik.full + log(N.ex)*npar.full + npar.full
  fullvec <- c(loglik.full, npar.full, AIC.full, BIC.full, cAIC.full)

  #------------ MML -----------
  N <- dim(object$X)[1]
  rv <- rowSums(object$X, na.rm = TRUE)                       #person raw scores
  npar.mml <- (dim(object$W)[2])#+(length(table(rv)))
  lmml <- sum(table(rv)*log(table(rv)/N))+object$loglik.cml   #MML likelihood
  AIC.mml <- -2*lmml + 2*npar.mml
  BIC.mml <- -2*lmml + log(N)*npar.mml
  cAIC.mml <- -2*lmml + log(N)*npar.mml + npar.mml
  mmlvec <- c(lmml, npar.mml, AIC.mml, BIC.mml, cAIC.mml)
  
  #------------- CML ---------------
  npar.cml <- dim(object$W)[2]
  lcml <- object$loglik.cml
  AIC.cml <- -2*lcml + 2*npar.cml
  BIC.cml <- -2*lcml + log(N)*npar.cml
  cAIC.cml <- -2*lcml + log(N)*npar.cml + npar.cml
  cmlvec <- c(lcml, npar.cml, AIC.cml, BIC.cml, cAIC.cml)
  
  ICtable <- rbind(fullvec, mmlvec, cmlvec)
  rownames(ICtable) <- c("joint log-lik", "marginal log-lik", "conditional log-lik")
  colnames(ICtable) <- c("value", "npar", "AIC", "BIC", "cAIC")

  result <- list(ICtable = ICtable)
  class(result) <- "ICr"
  result
}
