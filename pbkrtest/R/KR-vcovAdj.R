
## --------------------------------------------------------------------
## Calculate the adjusted covariance matrix for a mixed model
##
## Implemented in Banff, august 2013; Søren Højsgaard
## --------------------------------------------------------------------

vcovAdj <- function(object, details=0){
  UseMethod("vcovAdj")
}

vcovAdj.lmerMod <-
    vcovAdj.mer <-
        function(object, details=0){
    if (!(getME(object, "is_REML"))) {
        object <- update(object, . ~ ., REML = TRUE)
    }
    Phi      <- vcov(object)
    SigmaG   <- get_SigmaG( object, details )
    X        <- getME(object,"X")
    vcovAdj16_internal( Phi, SigmaG, X, details=details)
}


.vcovAdj_internal <- function(Phi, SigmaG, X, details=0){

    ##cat("vcovAdj_internal\n")
    ##SG<<-SigmaG
    DB <- details > 0 ## debugging only

    #print("HHHHHHHHHHHHHHH")
    #print(system.time({chol( forceSymmetric(SigmaG$Sigma) )}))
    #print(system.time({chol2inv( chol( forceSymmetric(SigmaG$Sigma) ) )}))

    ## print("HHHHHHHHHHHHHHH")
    ## Sig <- forceSymmetric( SigmaG$Sigma )
    ## print("HHHHHHHHHHHHHHH")
    ## print(system.time({Sig.chol <- chol( Sig )}))
    ## print(system.time({chol2inv( Sig.chol )}))

    t0 <- proc.time()
    ## print("HHHHHHHHHHHHHHH")
    SigmaInv <- chol2inv( chol( forceSymmetric(SigmaG$Sigma) ) )
    ## print("DONE --- HHHHHHHHHHHHHHH")

    if(DB){
        cat(sprintf("Finding SigmaInv: %10.5f\n", (proc.time()-t0)[1] ));
        t0 <- proc.time()
    }
    ##print("iiiiiiiiiiiii")

    t0 <- proc.time()
    ## Finding, TT, HH, 00
    n.ggamma <- SigmaG$n.ggamma
    TT       <- SigmaInv %*% X
    HH       <- OO <- vector("list", n.ggamma)
    for (ii in 1:n.ggamma) {
        .tmp <- SigmaG$G[[ii]] %*% SigmaInv
        HH[[ ii ]] <- .tmp
        OO[[ ii ]] <- .tmp %*% X
    }
    if(DB){cat(sprintf("Finding TT,HH,OO  %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
    ## if(DB){
    ##     cat("HH:\n"); print(HH); HH <<- HH
    ##     cat("OO:\n"); print(OO); OO <<- OO
    ## }

    ## Finding PP, QQ
    PP <- QQ <- NULL
    for (rr in 1:n.ggamma) {
        OrTrans <- t( OO[[ rr ]] )
        PP <- c(PP, list(forceSymmetric( -1 * OrTrans %*%  TT)))
        for (ss in rr:n.ggamma) {
            QQ <- c(QQ,list(OrTrans %*% SigmaInv %*% OO[[ss]] ))
        }}
    if(DB){cat(sprintf("Finding PP,QQ:    %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
    ## if(DB){
    ##     cat("PP:\n"); print(PP); PP2 <<- PP
    ##     cat("QP:\n"); print(QQ); QQ2 <<- QQ
    ## }

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

    WW <- if(condi>1e-10) forceSymmetric(2* solve(IE2)) else forceSymmetric(2* ginv(IE2))

    ## print("vcovAdj")
    UU <- matrix(0, nrow=ncol(X), ncol=ncol(X))
    ## print(UU)
    for (ii in 1:(n.ggamma-1)) {
        for (jj in c((ii+1):n.ggamma)) {
            www <- .indexSymmat2vec( ii, jj, n.ggamma )
            UU <- UU + WW[ii,jj] * (QQ[[ www ]] - PP[[ii]] %*% Phi %*% PP[[jj]])
        }}
    ## print(UU)

    UU <- UU + t(UU)
    ## UU <<- UU
    for (ii in 1:n.ggamma) {
        www <- .indexSymmat2vec( ii, ii, n.ggamma )
        UU<- UU +   WW[ii,ii] * (QQ[[ www ]] - PP[[ii]] %*% Phi %*% PP[[ii]])
    }
    ## print(UU)
    GGAMMA <-  Phi %*% UU %*% Phi
    PhiA   <-  Phi + 2 * GGAMMA
    attr(PhiA, "P")     <-PP
    attr(PhiA, "W")     <-WW
    attr(PhiA, "condi") <- condi
    PhiA
}
