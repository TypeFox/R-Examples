
## --------------------------------------------------------------------
## Calculate the adjusted covariance matrix for a mixed model
## --------------------------------------------------------------------

vcovAdj0 <- function(object, details=0) {
  DB <- details > 0 ## debugging only

  if(!.is.lmm(object)) {
    cat("Error in vccovAdj\n")
    cat(sprintf("model is not a linear mixed moxed model fitted with lmer\n"))
    stop()
  }

  if (!(getME(object, "is_REML"))){    #cat("\n object has been refitted with REML=TRUE \n")
    object <- update(object, .~., REML=TRUE)
  }

  ## Ready to go...

  X        <- getME(object,"X")
  Phi      <- vcov(object)
  SigmaG   <- LMM_Sigma_G( object, details )



  SigmaInv <- chol2inv( chol( forceSymmetric(SigmaG$Sigma) ) )
  if(DB){cat(sprintf("Finding SigmaInv: %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  ## Finding, TT, HH, 00
  n.ggamma <- SigmaG$n.ggamma
  TT       <- SigmaInv %*% X
  HH       <- OO <-  vector("list", n.ggamma)
  for (ii in 1:n.ggamma) {
    .DUM<-SigmaG$G[[ii]] %*% SigmaInv
    HH[[ii]] <- .DUM
    OO[[ii]] <- .DUM %*% X
  }
  if(DB){cat(sprintf("Finding TT,HH,OO  %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  ## Finding PP, QQ
  PP <- QQ <- NULL
  for (rr in 1:n.ggamma) {
    OrTrans <- t(OO[[rr]])
    PP <- c(PP, list(forceSymmetric( -1 * OrTrans %*%  TT)))
    for (ss in rr:n.ggamma) {
      QQ <- c(QQ,list(OrTrans %*% SigmaInv %*% OO[[ss]] ))
    }}
  if(DB){cat(sprintf("Finding PP,QQ:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  Ktrace <- matrix( NA, nrow=n.ggamma, ncol=n.ggamma )
  for (rr in 1:n.ggamma) {
    HrTrans <- t( HH[[rr]] )
    for (ss in rr:n.ggamma){
      Ktrace[rr,ss] <- Ktrace[ss,rr]<- sum( HrTrans * HH[[ss]] )
    }}
  if(DB){cat(sprintf("Finding Ktrace:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  ## Finding information matrix
  IE2 <- matrix( NA, nrow=n.ggamma, ncol=n.ggamma )
  for (ii in 1:n.ggamma) {
    Phi.P.ii <- Phi %*% PP[[ii]]
    for (jj in c(ii:n.ggamma)) {
      www <- .indexSymmat2vec( ii, jj, n.ggamma )

      IE2[ii,jj]<- IE2[jj,ii] <- Ktrace[ii,jj] -
        2 * sum(Phi*QQ[[ www ]]) + sum( Phi.P.ii * ( PP[[jj]] %*% Phi))
    }}
  if(DB){cat(sprintf("Finding IE2:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  eigenIE2 <- eigen(IE2,only.values=TRUE)$values
  condi    <- min(abs(eigenIE2))

  W <- if(condi>1e-10) forceSymmetric(2* solve(IE2)) else forceSymmetric(2* ginv(IE2))

  U <- matrix(0, nrow=ncol(X), ncol=ncol(X))
  for (ii in 1:(n.ggamma-1)) {
    for (jj in c((ii+1):n.ggamma)) {
      www <- .indexSymmat2vec( ii, jj, n.ggamma )
      U <- U + W[ii,jj] * (QQ[[ www ]] - PP[[ii]] %*% Phi %*% PP[[jj]])
    }}

### FIXME: Ulrich: Er det ikke sådan, at du her får beregnet diagonalen med to gange???
  U <- U + t(U)

  for (ii in 1:n.ggamma) {
    www <- .indexSymmat2vec( ii, ii, n.ggamma )
    U<- U +   W[ii,ii] * (QQ[[ www ]] - PP[[ii]] %*% Phi %*% PP[[ii]])
  }

  GGAMMA <-  Phi %*% U %*% Phi
  PhiA   <-  Phi + 2 * GGAMMA
  attr(PhiA, "P")     <-PP
  attr(PhiA, "W")     <-W
  attr(PhiA, "condi") <- condi
  PhiA
}



 ## Nov. 24. 2011; SHD
  ## Alternative computation of Ktrace. Seems to be no faster than the one above but please do
  ## not delete
  ##     Ktrace2 <- matrix(NA,n.ggamma,n.ggamma)
  ##     for (rr in 1:n.ggamma) {
  ##       HrTrans<-t(H[[rr]])
  ##       Ktrace2[rr,rr] <- sum( HrTrans * t(HrTrans))
  ##       if (rr < n.ggamma){
  ##         for (ss in (rr+1):n.ggamma) {
  ##           Ktrace2[rr,ss] <- Ktrace2[ss,rr]<- sum( HrTrans * H[[ss]])
  ##         }}}
  ##     cat(sprintf("Finding Ktrace(2): %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()
  ##     print(Ktrace)
  ##     print(Ktrace2)
  ##     print(Ktrace2-Ktrace)

