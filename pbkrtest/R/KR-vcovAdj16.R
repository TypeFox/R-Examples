.vcovAdj16 <-
    function(object, details=0){
        if (!(getME(object, "is_REML"))) {
            object <- update(object, . ~ ., REML = TRUE)
        }
        Phi      <- vcov(object)
        SigmaG   <- get_SigmaG( object, details )
        X        <- getME(object,"X")
        vcovAdj16_internal( Phi, SigmaG, X, details=details)
    }

vcovAdj16_internal <- function(Phi, SigmaG, X, details=0){

    details=0
    DB <- details > 0 ## debugging only
    t0 <- proc.time()
    
    ##Sigma <- SigmaG$Sigma
    n.ggamma <- SigmaG$n.ggamma

    M <- cbind(do.call(cbind, SigmaG$G), X)
    if(DB){cat(sprintf("M etc:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
    ##SinvM <- solve(SigmaG$Sigma, M, sparse=TRUE)
    SinvM <- chol2inv(chol( forceSymmetric( SigmaG$Sigma ))) %*% M
    ##SigmaInv <- chol2inv( chol( forceSymmetric(SigmaG$Sigma) ) )

    if(DB){cat(sprintf("SinvM etc:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
    
    v   <- c(rep(1:length(SigmaG$G), each=nrow(SinvM)), rep(length(SigmaG$G)+1, ncol(X)))
    idx <- lapply(unique.default(v), function(i) which(v==i))

    SinvG <- lapply(idx, function(z) SinvM[,z])  ## List of SinvG1, SinvG2,... SinvGr, SinvX
    SinvX <- SinvG[[length(SinvG)]]              ## Kaldes TT andre steder
    SinvG[length(SinvG)] <- NULL                 ## Er HH^t 

    if(DB){cat(sprintf("SinvG etc:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

    ##stat <<- list(SigmaG=SigmaG, X=X, M=M)
    
    OO <- lapply(1:n.ggamma, function(i) {
        SigmaG$G[[i]] %*% SinvX  ## G_i \Sigma\inv X; n \times p 
    })

    if(DB){cat(sprintf("Finding OO:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
    
    PP <- vector("list", n.ggamma)
    QQ <- vector("list", n.ggamma * (n.ggamma + 1) / 2 )
    index <- 1
    for (r in 1:n.ggamma) {
        OOt.r <- t( OO[[ r ]] )
        #str(list("dim(OOt.r)"=dim(OOt.r), "dim(SinvX)"=dim(SinvX)))
        ##PP[[r]] <- forceSymmetric( -1 * OOt.r %*%  SinvX) ## PP : p \times p
        PP[[r]] <- -1 * (OOt.r %*%  SinvX) ## PP : p \times p
        
        for (s in r:n.ggamma) {
            QQ[[index]] <- OOt.r %*% ( SinvG[[s]] %*% SinvX )
            index <- index + 1;
        }
    }
    ##stat16 <<- list(Phi=Phi, OO=OO, PP=PP,QQ=QQ)

    if(DB){cat(sprintf("Finding PP,QQ:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
    
    
    Ktrace <- matrix( NA, nrow=n.ggamma, ncol=n.ggamma )
    for (r in 1:n.ggamma) {
        HHr <- SinvG[[r]]
        for (s in r:n.ggamma){
            Ktrace[r,s] <- Ktrace[s,r] <- sum( HHr * SinvG[[s]] )
        }}

    if(DB){cat(sprintf("Finding Ktrace:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
        
    ## Finding information matrix
    IE2 <- matrix(0, nrow=n.ggamma, ncol=n.ggamma )
    for (ii in 1:n.ggamma) {
        Phi.P.ii <- Phi %*% PP[[ii]]
        for (jj in c(ii:n.ggamma)) {
            www <- .indexSymmat2vec( ii, jj, n.ggamma )
            IE2[ii,jj]<- IE2[jj,ii] <- Ktrace[ii,jj] -
                2 * sum(Phi * QQ[[ www ]]) + sum( Phi.P.ii * ( PP[[jj]] %*% Phi))
        }}
    if(DB){cat(sprintf("Finding IE2:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

    eigenIE2 <- eigen( IE2, only.values=TRUE )$values
    condi    <- min( abs( eigenIE2 ) )
    WW <- if ( condi > 1e-10 )
              forceSymmetric(2 * solve(IE2))
          else
              forceSymmetric(2 * ginv(IE2))
    
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
    for (ii in 1:n.ggamma) {
        www <- .indexSymmat2vec( ii, ii, n.ggamma )
        UU  <- UU + WW[ii,ii] * (QQ[[ www ]] - PP[[ii]] %*% Phi %*% PP[[ii]])
    }
    if(DB){cat(sprintf("Finding UU:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}    
    ## print(UU)
    GGAMMA <-  Phi %*% UU %*% Phi
    PhiA   <-  Phi + 2 * GGAMMA
    attr(PhiA, "P")     <- PP
    attr(PhiA, "W")     <- WW
    attr(PhiA, "condi") <- condi
    PhiA
}

