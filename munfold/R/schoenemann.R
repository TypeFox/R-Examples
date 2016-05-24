#Metric Unfolding after:
#Schoenemann, Peter E. 1970.
#On Metric Multidimensional Unfolding.
#Psychometrika 35,5:349-366
#
#Author: Martin Elff
#(c) 2002
#
#
#


schoenemann.unfolding <- function(   D,
                                        m=NULL,
                                        svd.eps=1e-7,
                                        eigen.eps=1e-7
                                        # w=NULL,
                                    ){
# Schoenemann's Step [1]: Read in data
    D <- as.matrix(D)
    if(ncol(D)>nrow(D)) {
            D<-t(D)
            D.transposed <- TRUE
        }
    else D.transposed <- FALSE
    if(min(D)<0) warning("Negative distances in argument matrix")
    p <- nrow(D)
    q <- ncol(D)
# Step [2]: Compute "quasi-scalar products matrix"
    D.rowmean <- apply(D,1,mean)
    D.colmean <- apply(D,2,mean)
    D.totmean <- mean(D)
    C <- sweep(D,1,D.rowmean)
    C <- sweep(C,2,D.colmean)
    C <- C + D.totmean
    C <- -0.5*C
# Step [3]: Decide on dimensionality and factor C
    S <- svd(C) # SVD is an Eckart-Young decomposition ...
    if(is.null(m))
        m <- length(S$d[S$d >= svd.eps])
    if(length(S$d[S$d >= svd.eps]) < m){
        if (length(S$d[S$d > 0]) <m) {
            warning("Dimension m of unfolding space reduced to the number of positive singular values of C")
            m <- length(S$d[S$d > 0])
        }
        else {
            m <- length(S$d[S$d >= svd.eps])
            warning("The number of singular values of C large enough is only ",m)
        }
    }
    if((p-1) < m*(m+3)/2)
        stop("Problem not identified - please specify a lower dimension (m)!")
    G <-S$u[,1:m] %*% diag(sqrt(S$d[1:m]),nrow=length(S$d[1:m]))
    H <-S$v[,1:m] %*% diag(sqrt(S$d[1:m]),nrow=length(S$d[1:m]))
# Step [4]: Construct coefficient matrix K_0
    K <- matrix(0,nrow=p,ncol=0)
    for(r in 1:m)
        for(s in r:m)
            K <- cbind(K, G[,r] * G[,s])
    for(r in 1:m)
        K <- cbind(K, G[,r])
    K <- sweep(K[1:(p-1),],2,K[p,])
    K[] <- ifelse(abs(K)<=100*.Machine$double.eps,0,K)
# Step [5]: Construct the vector of constants \bar{f}
    Fm <- D + 2*C
    Fm <- sweep(Fm[1:(p-1),],2,Fm[p,])
    Fm[] <- ifelse(abs(Fm)<=100*.Machine$double.eps,0,Fm)
    fb <- apply(Fm,1,mean)
    if(all(fb==0)) fb.null <- TRUE
    else fb.null <- FALSE
    if(!fb.null) {
    # Step [6]: Solve the K%*%unknown = fbar in least square sense
        fit <- lm.fit(x=K,y=fb)
        # Sometimes K is not full collumn rank ...
        if(fit$rank<NCOL(K)) {
        #    fit <- lm.fit(x=K[,1:fit$rank],y=fb)
            K.full.rank <- FALSE
            warning("K is not full collumn rank")
            }
        else K.full.rank <- TRUE
        nu <- fit$coefficients[1:((m*(m+1))/2)]
        nu <- ifelse(is.finite(nu),nu,0)
        xi <- fit$coefficients[(1+(m*(m+1))/2):((m*(m+3))/2)]
        xi <- ifelse(is.finite(xi),xi,0)
    # Step [7]: Construct matrix M from nu
        M <- matrix(0,nrow=m,ncol=m)
        ii <- 1
        for(r in 1:m)
            for(s in r:m) {
                if(r==s) M[r,s] <- nu[ii]
                else M[r,s] <- M[s,r] <- nu[ii]/2
                ii <- ii + 1
            }
    # Step [8]: Check definiteness of M
        M.eigen <- eigen(M)
        M.eigenvalues <-M.eigen$values
        if(all(M.eigen$values<=0)){
                warning("Only non-positive eigenvalues present in M, reflecting them")
                M.eigen$values <- -M.eigen$values
                eigenvalues.ok <- FALSE
        }
        else eigenvalues.ok <- TRUE
        if(min(M.eigen$values)<=0){
                warning(length(M.eigen$values<=0), " non-positive eigenvalues of M replaced by eigen.eps=",eigen.eps)
                M.eigen$values[M.eigen$values<=0] <- eigen.eps
        }
    # Step [9]: Factor M into TT'
        sqrt.eigen.values <- sqrt(M.eigen$values)
            Tm <- M.eigen$vectors %*% diag(sqrt.eigen.values,nrow=length(sqrt.eigen.values))
            Tm.inv <- diag(1/sqrt.eigen.values,nrow=length(sqrt.eigen.values)) %*% t(M.eigen$vectors)
    # Step [10]: Recover A_0
        A0 <- G %*% Tm
    # Step [11]: Recover the difference a0 - b0
        diffa0b0 <- Tm.inv %*% (xi/2)
    # Step [12]: Recover B_0
        B0 <- sweep(H %*% t(Tm.inv),2,diffa0b0)
    # Finished!

    ## kludge!
    if(!K.full.rank || !eigenvalues.ok){
        warning("Result shows only within-group configs")
        A0 <- G #%*% diag(S$d[1:m]),nrow=length(S$d[1:m]))
        B0 <- H #%*% diag(S$d[1:m],nrow=length(S$d[1:m]))
        }
    }
    else {
        warning("Difference between D and -2C vanishes, skipping steps 6-12")
        A0 <- G #%*% diag(S$d[1:m]),nrow=length(S$d[1:m]))
        B0 <- H #%*% diag(S$d[1:m],nrow=length(S$d[1:m]))
    }
    if(!is.null(colnames(D))) rownames(B0) <- colnames(D)
    if(!is.null(rownames(D))) rownames(A0) <- rownames(D)
    Delta <- matrix(0,nrow=nrow(D),ncol=ncol(D))
        for(i in 1:m)
        Delta <- Delta + outer(A0[,i],B0[,i],"-")^2
    stress <- sum((sqrt(D)-sqrt(Delta))^2)/sum(D)

    return(
        structure(
            list(
                A=A0,
                B=B0,
                d=S$d,
                fitted=Delta,
                stress=stress,
                dim=m
            ),
            class="unfolding"
        )
    )
}


