`probtrans` <- function(object,predt,direction=c("forward","fixedhorizon"),
                method=c("aalen","greenwood"),variance=TRUE,covariance=FALSE)
{
    if (!inherits(object, "msfit"))
        stop("'object' must be a 'msfit' object")
    method <- match.arg(method)
    direction <- match.arg(direction)
    trans <- object$trans
    transit <- to.trans2(trans)
    numtrans <- nrow(transit)
    stackhaz <- object$Haz
    stackvarhaz <- object$varHaz
    for (i in 1:numtrans)
        stackhaz$dhaz[stackhaz$trans==i] <- diff(c(0,stackhaz$Haz[stackhaz$trans==i]))
    if (direction=="forward")
        stackhaz <- stackhaz[stackhaz$time > predt,]
    else stackhaz <- stackhaz[stackhaz$time <= predt,]

    untimes <- sort(unique(stackhaz$time))
    TT <- length(untimes)
    S <- nrow(trans)

    if (covariance) variance <- TRUE # if covariance=TRUE and variance=FALSE, variance is overruled
    if (direction=="forward") {
        if (variance==TRUE) res <- array(0,c(TT+1,2*S+1,S)) # 2*S+1 for time, probs (S), se (S)
        else res <- array(0,c(TT+1,S+1,S)) # S+1 for time, probs (S)
        # first line (in case of forward) contains values for t=predt
        res[1,1,] <- predt
        for (j in 1:S) res[1, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
        if (variance) res[1,(S+2):(2*S+1),] <- 0
    }
    else {
        # situation for backward is different from forward,
        # depends on whether predt is an event time
        if (predt %in% untimes) {
            if (variance) res <- array(0,c(TT+1,2*S+1,S)) # t=0 and all event times
            else res <- array(0,c(TT+1,S+1,S))
            res[TT+1,1,] <- predt
            for (j in 1:S) res[TT+1, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
            if (variance) res[TT+1,(S+2):(2*S+1),] <- 0
        }
        else {
            if (variance) res <- array(0,c(TT+2,2*S+1,S)) # t=0, event times, t=predt
            else res <- array(0,c(TT+2,S+1,S))
            res[TT+1,1,] <- max(untimes)
            for (j in 1:S) res[TT+1, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
            if (variance) res[TT+1,(S+2):(2*S+1),] <- 0
            res[TT+2,1,] <- predt
            for (j in 1:S) res[TT+2, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
            if (variance) res[TT+2,(S+2):(2*S+1),] <- 0
        }
    }
    P <- diag(S)
    
    if (covariance) {
      varParr <- array(0,c(S^2,S^2,TT+1))
      if ((direction=="fixedhorizon") & !(predt %in% untimes))
        varParr <- array(0,c(S^2,S^2,TT+2))
      ffrom <- rep(1:S, S)
      tto <- rep(1:S, rep(S,S))
      fromto <- paste("from",ffrom,"to",tto,sep="")
      if (direction=="forward")
        dimnames(varParr) <- list(fromto,fromto,c(predt,untimes))
      else {
        if (predt %in% untimes) dimnames(varParr) <- list(fromto,fromto,c(0,untimes))
        else dimnames(varParr) <- list(fromto,fromto,c(0,untimes,predt))
      }
    }
    if (variance) {
        varP <- matrix(0,S^2,S^2)
        if (direction=="forward") {
            varAnew <- array(0,c(S,S,S,S))
            if (predt !=0) {
                tmin <- max(stackvarhaz$time[stackvarhaz$time <=predt])
                varHaz <- stackvarhaz[stackvarhaz$time==tmin,]
                lHaz <- nrow(varHaz)
                for (j in 1:lHaz)
                {
                    from1 <- transit$from[transit$transno==varHaz$trans1[j]]
                    to1 <- transit$to[transit$transno==varHaz$trans1[j]]
                    from2 <- transit$from[transit$transno==varHaz$trans2[j]]
                    to2 <- transit$to[transit$transno==varHaz$trans2[j]]
                    varAnew[from1, to1, from2, to2] <-
                        varAnew[from2, to2, from1, to1] <- varHaz$varHaz[j]
                }
            }
        }
        else { 
            # varA on last event time (only those elements borrowed
            # directly from object, rest is calculated later)
            varA <- array(0,c(S,S,S,S))
            varHaz <- stackvarhaz[stackvarhaz$time==untimes[TT],]
            lHaz <- nrow(varHaz)
            for (j in 1:lHaz)
            {
                from1 <- transit$from[transit$transno==varHaz$trans1[j]]
                to1 <- transit$to[transit$transno==varHaz$trans1[j]]
                from2 <- transit$from[transit$transno==varHaz$trans2[j]]
                to2 <- transit$to[transit$transno==varHaz$trans2[j]]
                varA[from1, to1, from2, to2] <-
                    varA[from2, to2, from1, to1] <- varHaz$varHaz[j]
            }
        }
    }
    
    for (i in 1:TT)
    {
        idx <- ifelse(direction=="forward",i,TT+1-i)
        tt <- untimes[idx]
        Haztt <- stackhaz[stackhaz$time==tt,]
        lHaztt <- nrow(Haztt)
        # build S x S matrix IplusdA
        IplusdA <- diag(S)
        for (j in 1:lHaztt)
        {
            from <- transit$from[transit$transno==Haztt$trans[j]]
            to <- transit$to[transit$transno==Haztt$trans[j]]
            IplusdA[from, to] <- Haztt$dhaz[j]
            IplusdA[from, from] <- IplusdA[from, from] - Haztt$dhaz[j]
        }
        if (any(diag(IplusdA)<0))
            warning("Warning! Negative diagonal elements of (I+dA); the estimate may not be meaningful. \n")
        
        if (variance) {
            if (direction=="forward"){
                varA <- varAnew
                varAnew <- array(0,c(S,S,S,S))
                varHaztt <- stackvarhaz[stackvarhaz$time==tt,]
                lHaztt <- nrow(varHaztt)
                for (j in 1:lHaztt)
                {
                    from1 <- transit$from[transit$transno==varHaztt$trans1[j]]
                    to1 <- transit$to[transit$transno==varHaztt$trans1[j]]
                    from2 <- transit$from[transit$transno==varHaztt$trans2[j]]
                    to2 <- transit$to[transit$transno==varHaztt$trans2[j]]
                    varAnew[from1, to1, from2, to2] <-
                        varAnew[from2, to2, from1, to1] <- varHaztt$varHaz[j]
                }
                vardA <- varAnew - varA
            }
            else {
                varAttmin <- array(0,c(S,S,S,S))
                varHazttmin <- stackvarhaz[stackvarhaz$time==untimes[idx-1],]
                lHazttmin <- nrow(varHazttmin)
                for (j in 1:lHazttmin)
                {
                    from1 <- transit$from[transit$transno==varHazttmin$trans1[j]]
                    to1 <- transit$to[transit$transno==varHazttmin$trans1[j]]
                    from2 <- transit$from[transit$transno==varHazttmin$trans2[j]]
                    to2 <- transit$to[transit$transno==varHazttmin$trans2[j]]
                    varAttmin[from1, to1, from2, to2] <-
                        varAttmin[from2, to2, from1, to1] <- varHazttmin$varHaz[j]
                }
                vardA <- varA - varAttmin
                varA <- varAttmin # ready for the next round
            }
        
            for (from in 1:S)
            {
                for (from2 in 1:S)
                {
                    for (to2 in 1:S)
                    {
                        if (to2!=from2)
                            vardA[from, from, from2, to2] <-
                                vardA[from2, to2, from, from] <-
                                    -sum(vardA[from,-from,from2,to2])
                    }
                }
            }
            for (from in 1:S)
            {
                for (from2 in 1:S)
                    vardA[from,from,from2,from2] <-
                        vardA[from2,from2,from,from] <-
                        -sum(vardA[from,from,from2,-from2])
            }
            vardA <- matrix(vardA,S^2,S^2)
        }
        
        if (method=="aalen")
        {
            if (direction=="forward") 
            {  
                P <- P %*% IplusdA
                if (variance) {
                    tmp1 <- kronecker(t(IplusdA),diag(S)) %*% varP %*% kronecker(IplusdA,diag(S))
                    tmp2 <- kronecker(diag(S),P) %*% vardA %*% kronecker(diag(S),t(P))
                    varP <- tmp1 + tmp2
                }
            }
            else {
                if (variance) {
                    tmp1 <- kronecker(diag(S), IplusdA) %*% varP %*% kronecker(diag(S), t(IplusdA))
                    tmp2 <- kronecker(t(P), IplusdA) %*% vardA %*% kronecker(P,t(IplusdA))
                    varP <- tmp1+tmp2
                }
                P <- IplusdA %*% P
            }
        }
        if (method=="greenwood")
        {
            if (direction=="forward")
            {
                if (variance) {
                    tmp1 <- kronecker(t(IplusdA),diag(S)) %*% varP %*% kronecker(IplusdA,diag(S))
                    tmp2 <- kronecker(diag(S),P) %*% vardA %*% kronecker(diag(S),t(P))
                    varP <- tmp1 + tmp2
                }
                P <- P %*% IplusdA
            }
            else {
                if (variance) {
                    tmp1 <- kronecker(diag(S), IplusdA) %*% varP %*% kronecker(diag(S), t(IplusdA))
                    tmp2 <- kronecker(t(P), diag(S)) %*% vardA %*% kronecker(P, diag(S))
                    varP <- tmp1+tmp2
                }
                P <- IplusdA %*% P
            }
        }
        if (variance) {
            seP <- sqrt(diag(varP))
            seP <- matrix(seP,S,S)
        }
        if (covariance) {
          if (direction=="forward") varParr[,,i+1] <- varP
          else {
            varParr[,,idx+1] <- varP
          }
        }
        if (direction=="forward") {
            res[idx+1,1,] <- tt
            res[idx+1,2:(S+1),] <- t(P)
            if (variance) res[idx+1,(S+2):(2*S+1),] <- t(seP)
        }
        else {
            res[idx,1,] <- ifelse(i==TT,0,untimes[TT-i])
            res[idx,2:(S+1),] <- t(P)
            if (variance) res[idx,(S+2):(2*S+1),] <- t(seP)
        }
    }
    if (covariance & (direction=="fixedhorizon"))
      varParr[,,1] <- varParr[,,2]
    ### res[,,s] contains prediction from state s, convert to list of dataframes
    res2  <- vector("list", S)
    for (s in 1:S) {
        tmp <- as.data.frame(res[,,s])
        if (min(dim(tmp))==1) tmp <- res[,,s]
        if (variance) names(tmp) <- c("time",paste("pstate",1:S,sep=""),paste("se",1:S,sep=""))
        else names(tmp) <- c("time",paste("pstate",1:S,sep=""))
        res2[[s]] <- tmp
    }
    if (covariance) res2$varMatrix <- varParr
    res2$trans <- trans
    res2$method <- method
    res2$predt <- predt
    res2$direction <- direction
    class(res2) <- "probtrans"
    return(res2)
}
