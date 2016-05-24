FLLat.FDR <- function(Y,Y.FLLat,n.thresh=50,fdr.control=0.05,pi0=1,
                      n.perms=20) {

    ## Error checking.
    if (!all(is.matrix(Y),is.double(Y))) {
        stop("'Y' must be a numeric matrix")
    }
    if (!inherits(Y.FLLat,"FLLat")) {
        stop("'Y.FLLat' must be of class 'FLLat'")
    }
    if (!all(length(n.thresh)==1,n.thresh>0)) {
        stop("'n.thresh' must be an integer > 0")
    }
    if (!all(length(fdr.control)==1,fdr.control>0,fdr.control<=1)) {
        stop("'fdr.control' must be a real number between 0 and 1")
    }
    if (!all(length(pi0)==1,pi0>0,pi0<=1)) {
        stop("'pi0' must be a real number between 0 and 1")
    }
    if (!all(length(n.perms)==1,n.perms>0)) {
        stop("'n.perms' must be an integer > 0")
    }

    ## S, L and J.
    S <- ncol(Y)
    L <- nrow(Y)
    J <- ncol(Y.FLLat$Beta)

    ## Lambda1 and lambda2.
    lam1 <- Y.FLLat$lam1; lam2 <- Y.FLLat$lam2

    ## Threshold values for which to estimate FDRs.
    pos.fit <- abs(Y.FLLat$Beta%*%Y.FLLat$Theta) 
    thresh.vals <- seq(0,max(pos.fit),len=n.thresh+1)[-1]

    ## The denominator values for the FDR estimates.
    fdr.dens <- rep(NA,length(thresh.vals))
    for (i in 1:length(thresh.vals)) {
        fdr.dens[i] <- sum(pos.fit>=thresh.vals[i])
    }

    ## Matrix of FDRs for the permutations.
    fdrs <- matrix(0,nrow=length(thresh.vals),ncol=n.perms)

    for (i in 1:n.perms) {

        ## The permuted data.
        perm.Y <- matrix(NA,nrow=L,ncol=S)
        for (j in 1:S) {
            perm.Y[,j] <- Y[sample(L),j]
        }

        ## FLLat.
        perm.FLLat <- FLLat(Y=perm.Y,J=J,lam1=lam1,lam2=lam2)

        ## Positive fitted values.
        perm.pos.fit <- abs(perm.FLLat$Beta%*%perm.FLLat$Theta)

        ## FDRs.
        for (j in 1:length(thresh.vals)) {
            if (fdr.dens[j]==0) {
                fdrs[j,i] <- 0
            } else {
                fdrs[j,i] <- sum(perm.pos.fit>=thresh.vals[j])/fdr.dens[j]
            }
        }
    }

    ## Average FDRs.
    avg.fdrs <- rowMeans(fdrs)*pi0

    ## Threshold which controls FDR at fdr.control.
    thresh.control <- min(thresh.vals[avg.fdrs<=fdr.control])

    results <- list("thresh.values"=thresh.vals,"FDRs"=avg.fdrs,
                    "thresh.control"=thresh.control)
    class(results) <- "FDR"
    
    return(results)

}
