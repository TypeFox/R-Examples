compromise <- function(nfactors, G1, class=3, msg=TRUE){
    if (!class %in% c(1,2,3,4)) stop("class must be an integer from 1 to 4")
    if (!is.numeric(G1)) stop("G1 must be numeric")
    if (length(G1)==0) stop("At least one effect must be in group 1")
    if (length(G1)==nfactors) stop("Use a resolution V design")
    if (length(G1)==1 & class==1) stop("no estimable effects required!")
    if (!is.numeric(nfactors)) stop("nfactors must be numeric")
    if (!all(G1%%1==0)) stop("G1 must consist of integer values")
    if (!all(nfactors%%1==0)) stop("nfactors must be integer")
    if (!length(unique(G1))==length(G1)) stop("non-unique values in G1")
    if (!nfactors>=max(G1)) stop("G1 must not contain numbers larger than nfactors")
    if (!min(G1)>=1) stop("G1 must not contain numbers smaller than 1")
    G2 <- setdiff(1:nfactors,G1)
    if (length(G2)==0) stop("At least one effect must be in group 2")
    
    ## maximum number of factors for resolution V
    Mk <- c(2,3,5,6,8,11,17,23,32,47,65)
    k.cand <- 2:12
    names(Mk) <- k.cand   ## names are k (dimension of generating full factorial)
    
    perms <- combn(nfactors, length(G1))
    perms.full <- matrix(NA, ncol(perms), nfactors)
    for (i in 1:ncol(perms))
        perms.full[i,] <- c(perms[,i],setdiff(1:nfactors,perms[,i])) 
    perms.full <- perms.full[,invperm(c(G1,G2)),drop=FALSE]
    perms.full <- perms.full[ord(perms.full),,drop=FALSE]
    if (length(G1)>1)
    requirement <- apply(matrix(Letters[G1[combn(length(G1),2)]],nrow=2),2,"paste",collapse="")
    else requirement <- character(0)
    if (class==3)
    requirement <- c(requirement, outer(Letters[G1],Letters[G2],FUN=function(X,Y) paste(pmin(X,Y),pmax(X,Y),sep="")))
    if (class==2 & length(G2)>=2){
    requirement <- c(requirement, apply(matrix(Letters[G2[combn(length(G2),2)]],nrow=2),
                            2,"paste",collapse=""))
                            }
    if (class==4)
    requirement <- c(outer(Letters[G1],Letters[G2],FUN=function(X,Y) paste(pmin(X,Y),pmax(X,Y),sep="")))
    ## check existence rules for clear designs by Ke, Tang and Wu
    message <- NULL
    if (class==2) {
        minnrun.clear <- 2^k.cand[min(which(Mk>=nfactors))]
        message <- paste(gettext("a clear design requires at least"), 
              minnrun.clear,gettext("runs (resolution V)"))
          }
    else{
        check.0 <- 4*nfactors^2 + 4*nfactors - 2^(k.cand+3) + 9
        ## class 1
        if (class==1) {
          sel <- Mk-2 >= length(G1)
          if (any(check.0>=0))
               sel[check.0>=0] <- sel[check.0 >= 0] & 
                    nfactors - 1/2 - 1/2*sqrt(check.0[check.0>=0]) >= length(G1)
                    }
        ## class 3
        if (class==3){ 
             sel <- Mk-3 >= length(G1)
          if (any(check.0>=0))
               sel[check.0>=0] <- sel[check.0 >= 0] & 
                    nfactors - 3/2 - 1/2*sqrt(check.0[check.0>=0]) >= length(G1)
                    }
        ## class 4
        if (class==4){
             m1 <- min(length(G1), nfactors-length(G1))
             sel <- 2^k.cand >= m1
             check.0 <- nfactors^2 + 8*nfactors - 2^(k.cand+2) -4
             if (any(check.0>=0))
                 sel[check.0>=0] <- nfactors/2 - 1/2*sqrt(check.0[check.0>=0]) >= m1
        }
        pos <- min(c(which(sel), which(Mk >=nfactors)) )  ## Ke et al. rule or resolution V
        minnrun.clear <- 2^k.cand[pos]
        if (Mk[pos] >= nfactors) 
            message <- paste(gettext("a clear design requires at least"), minnrun.clear, gettext("runs (resolution V)"))
        else 
            message <- paste(gettext("a clear design requires at least"), minnrun.clear, gettext("runs"))
    }
    if (msg) message(message)
    list(perms.full=perms.full, requirement=requirement, class=class, minnrun.clear=minnrun.clear)
}
