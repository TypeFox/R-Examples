##setGeneric("PcaHubert", function(x, ...) standardGeneric("PcaHubert"))
##setMethod("PcaHubert", "formula", PcaHubert.formula)
##setMethod("PcaHubert", "ANY", PcaHubert.default)

setMethod("getQuan", "PcaHubert", function(obj) obj@quan)

##  The S3 version
PcaHubert <- function (x, ...) UseMethod("PcaHubert")

PcaHubert.formula <- function (formula, data = NULL, subset, na.action, ...)
{
    cl <- match.call()

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a `standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)

    res <- PcaHubert.default(x, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("PcaHubert")
    res@call <- cl

#    if (!is.null(na.act)) {
#        res$na.action <- na.act
#        if (!is.null(sc <- res$x))
#            res$x <- napredict(na.act, sc)
#    }

    res
}

PcaHubert.default <- function(x, k=0, kmax=10, alpha=0.75, mcd=TRUE, maxdir=250,
    scale=FALSE, signflip=TRUE, crit.pca.distances=0.975, trace=FALSE, ...)
{
## k    -   Number of principal components to compute. If \code{k} is missing,
##              or \code{k = 0}, the algorithm itself will determine the number of
##              components by finding such \code{k} that \code{l_k/l_1 >= 10.E-3} and
##              \code{sum_1:k l_j/sum_1:r l_j >= 0.8}. It is preferable to
##              investigate the scree plot in order to choose the number
##              of components and the run again. Default is \code{k=0}.
##
## kmax -   Maximal number of principal components to compute Default is \code{kmax=10}.
##              If \code{k} is provided, \code{kmax} does not need to be specified,
##              unless \code{k} is larger than 10.
## alpha    This parameter measures the fraction of outliers the algorithm should
##              resist. In MCD alpha controls the size of the subsets over which the determinant
##              is minimized, i.e. alpha*n observations are used for computing the determinant.
##              Allowed values are between 0.5 and 1 and the default is 0.5.
## mcd  -   Logical - when the number of variables is sufficiently small,
##              the loadings are computed as the eigenvectors of the MCD covariance matrix,
##              hence the function \code{\link{CovMcd}()} is automatically called. The number of
##              principal components is then taken as k = rank(x). Default is \code{mcd=TRUE}.
##              If \code{mcd=FALSE}, the ROBPCA algorithm is always applied.
## trace    whether to print intermediate results. Default is \code{trace = FALSE}}

##  Example:
##  data(hbk)
##  pca <- PcaHubert(hbk)

##
## The value returned by 'PcaHubert' is an S4 object containing the following slots:
##
## loadings     -   Robust loadings (eigenvectors)
## eigenvalues  -   Robust eigenvalues
## center       -   Robust center of the data
## scores       -   Robust scores
## k            -   Number of (chosen) principal components
##
## quan         -   The quantile h used throughout the algorithm
## sd           -   Robust score distances within the robust PCA subspace
## od           -   Orthogonal distances to the robust PCA subspace
## cutoff       -   Cutoff values for the robust score and orthogonal distances
## flag         -   The observations whose score distance is larger than result.cutoff.sd
##                  or whose orthogonal distance is larger than result$cutoff$od
##                  can be considered as outliers and receive a flag equal to zero.
##                  The regular observations receive a flag 1.

## This implementation followes closely the Matlab implementation, available as part of 'LIBRA,
##  a Matlab Library for Robust Analysis':
##      www.wis.kuleuven.ac.be/stat/robust.html

    cl <- match.call()

    if(missing(x))
        stop("You have to provide at least some data")

    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    ## ___Step 1___: Reduce the data space to the affine subspace spanned by the n observations
    ##  Apply svd() to the mean-centered data matrix. If n > p we use the kernel approach -
    ##  the decomposition is based on computing the eignevalues and eigenvectors of(X-m)(X-m)'

##    Xsvd <- kernelEVD(data, scale=scale, signflip=signflip)
    Xsvd <- classPC(data, scale=scale, signflip=signflip)

    if(Xsvd$rank == 0)
        stop("All data points collapse!")

    ## VT::27.08.2010: introduce 'scale' parameter; return the scale in the value object
    ##
    myscale = vector('numeric', p) + 1
    if(scale)
        myscale <- Xsvd$scale

    ##
    ## verify and set the input parameters: alpha, k and kmax
    ## determine h based on alpha and kmax, using the function h.alpha.n()
    ##

    ## VT::06.11.2012 - kmax <= floor(n/2) is too restrictive
    ##    kmax <- max(min(floor(kmax), floor(n/2), Xsvd$rank),1)
    kmax <- max(min(floor(kmax), Xsvd$rank),1)

    if((k <- floor(k)) < 0)
        k <- 0
    else if(k > kmax) {
        warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
        k <- kmax
    }

    if(missing(alpha))
    {
        default.alpha <- alpha
        h <- min(h.alpha.n(alpha, n, kmax), n)
        alpha <- h/n
        if(k == 0) {
            if(h < floor((n+kmax+1)/2)) {
                h <- floor((n+kmax+1)/2)
                alpha <- h/n
                warning(paste("h should be larger than (n+kmax+1)/2. It is set to its minimum value ", h, ".", sep=""))
            }
        }
        else {

            if(h < floor((n+k+1)/2)) {
                h <- floor((n+k+1)/2)
                alpha <- h/n
                warning(paste("h should be larger than (n+k+1)/2. It is set to its minimum value ", h, ".", sep=""))
            }
        }
        if(h > n) {

            alpha <- default.alpha
            if(k == 0)
                h <- h.alpha.n(alpha, n, kmax)
            else
                h <- h.alpha.n(alpha, n, k)
            warning(paste("h should be smaller than n = ", n, ". It is set to its default value ", h, ".", sep=""))
        }
    }else
    {

        if(alpha < 0.5 | alpha > 1)
            stop("Alpha is out of range: should be between 1/2 and 1")

        if(k == 0)
            h <- h.alpha.n(alpha, n, kmax)
        else
            h <- h.alpha.n(alpha, n, k)
    }

##    VT::30.10.2014 - classSVD and kernelEVD are simplified and
##      moved to robustbase, furthermore replaced by classPC:
##      no more return scores
##    X <- Xsvd$scores
    X <- scale(data, center=Xsvd$center, scale=Xsvd$scale) %*% Xsvd$loadings
    center <- Xsvd$center
    rot <- Xsvd$loadings

    ##
    ## ___Step 2___: Either calculate the standard PCA on the MCD covariance matrix (p<<n)
    ##  or apply the ROBPCA algorithm. If mcd=FALSE, allways apply ROBPCA.
    ##
    if(ncol(X) <= min(floor(n/5), kmax) & mcd)    # p << n => apply MCD
    {
        if(trace)
            cat("\nApplying MCD.\n")

        X.mcd <- CovMcd(as.data.frame(X), alpha=alpha)
        X.mcd.svd <- svd(getCov(X.mcd))

        rank <- ncol(X)     ## The covariance matrix is not singular
        ev <- X.mcd.svd$d

        ## VT::11.28.2015: Choose the number of components k (if not specified)
        ##
        ## Use the test l_k/l_1 >= 10.E-3, i.e. the ratio of
        ## the k-th eigenvalue to the first eigenvalue (sorted decreasingly) is larger than
        ## 10.E/3 and the fraction of the cumulative dispersion is larger or equal 80%
        ##
        if(k != 0)
            k <- min(k, ncol(X))
        else
        {
            test <- which(ev/ev[1] <= 1.E-3)
            k <- if(length(test) != 0)  min(min(rank, test[1]), kmax)
                 else                   min(rank, kmax)

            cumulative <- cumsum(ev[1:k])/sum(ev)
            if(cumulative[k] > 0.8) {
                k <- which(cumulative >= 0.8)[1]
            }

            if(trace)
                cat("\n k, kmax, rank, p: ", k, kmax, rank, p, "\n")
            if(trace)
                cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="")
        }

        scores <- (X - matrix(rep(getCenter(X.mcd), times=nrow(X)), nrow=nrow(X), byrow=TRUE)) %*% X.mcd.svd$u
        center <- as.vector(center + getCenter(X.mcd) %*% t(rot))
        eigenvalues <- X.mcd.svd$d[1:k]
        loadings <- Xsvd$loadings %*% X.mcd.svd$u[,1:k]
        scores <- as.matrix(scores[,1:k])
        if(is.list(dimnames(data)) && !is.null(dimnames(data)[[1]]))
        {
            dimnames(scores)[[1]] <- dimnames(data)[[1]]
        } else {
            dimnames(scores)[[1]] <- 1:n
        }
        dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))
        dimnames(scores)[[2]] <- as.list(paste("PC", seq_len(ncol(scores)), sep = ""))

        res <- new("PcaHubert",call=cl,
                            loadings=loadings,
                            eigenvalues=eigenvalues,
                            center=center,
                            scale=myscale,
                            scores=scores,
                            k=k,
                            quan=X.mcd@quan,
                            alpha=alpha,
                            n.obs=n)

    }
    else                                        # p > n or mcd=FALSE => apply the ROBPCA algorithm
    {
        if(trace)
            cat("\nApplying the projection method of Hubert.\n")

        ##
        ##  For each direction 'v' through 2 data points we project the n data points xi on v
        ##  and compute their robustly standardized absolute residual
        ##  |xi'v - tmcd(xj'v|/smcd(xj'v)
        ##
        alldir <- choose(n, 2)                  # all possible directions through two data points - n*(n-1)/2
        ndir <- min(maxdir, alldir)             # not more than maxdir (250)
        all <- (ndir == alldir)                 # all directions if n small enough (say n < 25)
        B <- extradir(X, ndir, all=all)         # the rows of B[ndir, ncol(X)] are the (random) directions

        Bnorm <- vector(mode="numeric", length=nrow(B))   # ndir x 1
        Bnorm <- apply(B, 1, vecnorm)           #
        Bnormr <- Bnorm[Bnorm > 1.E-12]         # choose only the nonzero length vectors
        m <- length(Bnormr)                     # the number of directions that will be used
        B <- B[Bnorm > 1.E-12,]                 # B[m x ncol(X)]
        A <- diag(1/Bnormr) %*% B               # A[m x ncol(X)]

        Y <- X %*% t(A)                         # n x m - projections of the n points on each of the m directions
        Z <- matrix(0,n, m)                     # n x m - to collect the outlyingness of each point on each direction
        for(i in 1:m) {
            umcd <- unimcd(Y[,i], quan=h)       # one-dimensional MCD: tmcd and smcd

            if(umcd$smcd < 1.E-12) {
                ## exact fit situation: will not be handled
                if((r2 <- rankMM(data[umcd$weights==1,])) == 1)
                    stop("At least ", sum(umcd$weights), " observations are identical.")
            }
            else
                Z[,i] <- abs(Y[,i] - umcd$tmcd) / umcd$smcd
        }

        H0 <- order(apply(Z, 1, max))           # n x 1 - the outlyingnesses of all n points
        Xh <- X[H0[1:h], ]                      # the h data points with smallest outlyingness

##        Xh.svd <- classSVD(Xh)
        Xh.svd <- classPC(Xh)
        kmax <- min(Xh.svd$rank, kmax)
        if(trace)
            cat("\nEigenvalues: ", Xh.svd$eigenvalues, "\n")

        ##
        ## Find the number of PC 'k'
        ## Use the test l_k/l_1 >= 10.E-3, i.e. the ratio of
        ## the k-th eigenvalue to the first eigenvalue (sorted decreasingly) is larger than
        ## 10.E/3 and the fraction of the cumulative dispersion is larger or equal 80%
        ##
        if(k == 0)
        {
            test <- which(Xh.svd$eigenvalues/Xh.svd$eigenvalues[1] <= 1.E-3)
            k <- if(length(test) != 0)  min(min(Xh.svd$rank, test[1]), kmax)
                 else                   min(Xh.svd$rank, kmax)

            cumulative <- cumsum(Xh.svd$eigenvalues[1:k])/sum(Xh.svd$eigenvalues)
            if(cumulative[k] > 0.8) {
                k <- which(cumulative >= 0.8)[1]
            }
            if(trace)
                cat(paste("The number of principal components is set by the algorithm. It is set to ", k, ".\n", sep=""))
        }

        if(trace)
            cat("\nXsvd$rank, Xh.svd$rank, k and kmax: ", Xsvd$rank, Xh.svd$rank, k, kmax,"\n")

        ## perform extra reweighting step
        if(k != Xsvd$rank)
        {
            ## VT::27.08.2010 - bug report from Stephen Milborrow: if n is small relative to p
            ## k can be < Xsvd$rank but larger than Xh.svd$rank - the number of observations in
            ## Xh.svd is roughly half of these in Xsvd
            k <- min(Xh.svd$rank, k)

##            XRc <- X - repmat(Xh.svd$center, n, 1)
            XRc <- X - matrix(rep(Xh.svd$center, times=n), nrow=n, byrow=TRUE)

            Xtilde <- XRc %*% Xh.svd$loadings[,1:k] %*% t(Xh.svd$loadings[,1:k])
            Rdiff <- XRc - Xtilde
            odh <- apply(Rdiff, 1, vecnorm)
            umcd <- unimcd(odh^(2/3), h)
            cutoffodh <- sqrt(qnorm(0.975, umcd$tmcd, umcd$smcd)^3)
            indexset <- (odh <= cutoffodh)

##            Xh.svd <- classSVD(X[indexset,])
            Xh.svd <- classPC(X[indexset,])
            k <- min(Xh.svd$rank, k)
        }

        ## Project the data points on the subspace spanned by the first k0 eigenvectors
        center <- center + Xh.svd$center %*% t(rot)
        rot <- rot %*% Xh.svd$loadings
#        X2 <- (X - repmat(Xh.svd$center, n, 1)) %*% Xh.svd$loadings
        X2 <- (X - matrix(rep(Xh.svd$center, times=n), nrow=n, byrow=TRUE)) %*% Xh.svd$loadings

        X2 <- as.matrix(X2[ ,1:k])
        rot <- as.matrix(rot[ ,1:k])

        ## Perform now MCD on X2 in order to obtain a robust scatter matrix: find h data points whose
        ##  covariance matrix has minimal determinant
        ##
        ## First apply C-step with the h points selected in the first step, i.e. those that
        ## determine the covariance matrix after the C-steps have converged.
        ##
        mah <- mahalanobis(X2, center=rep(0, ncol(X2)), cov=diag(Xh.svd$eigenvalues[1:k], nrow=k))
        oldobj <- prod(Xh.svd$eigenvalues[1:k])
        niter <- 100
        for(j in 1:niter) {
            if(trace)
                cat("\nIter=",j, " h=", h, " k=", k, " obj=", oldobj, "\n")
            Xh <- X2[order(mah)[1:h], ]

##            Xh.svd <- classSVD(as.matrix(Xh))
            Xh.svd <- classPC(as.matrix(Xh))
            obj <- prod(Xh.svd$eigenvalues)
#            X2 <- (X2 - repmat(Xh.svd$center, n, 1)) %*% Xh.svd$loadings
            X2 <- (X2 - matrix(rep(Xh.svd$center, times=n), nrow=n, byrow=TRUE)) %*% Xh.svd$loadings

            center <- center + Xh.svd$center %*% t(rot)
            rot <- rot %*% Xh.svd$loadings
            mah <- mahalanobis(X2, center=matrix(0,1, ncol(X2)), cov=diag(Xh.svd$eigenvalues, nrow=length(Xh.svd$eigenvalues)))
            if(Xh.svd$rank == k & abs(oldobj - obj) < 1.E-12)
                break

            oldobj <- obj
            if(Xh.svd$rank < k) {
                j <- 1
                k <- Xh.svd$rank
            }
        }

        ## Perform now MCD on X2
        X2mcd <- CovMcd(X2, nsamp=250, alpha=alpha)
        if(trace)
            cat("\nMCD crit=",X2mcd@crit," and C-Step obj function=",obj," Abs difference=", abs(X2mcd@crit-obj), "\n")

        ## VT::14.12.2009 - if there is even a slight difference between mcd$crit and obj
        ## and it is on the negative side, the following reweighting step will be triggered,
        ## which could lead to unwanted difference in the results. Therefore compare with
        ## a tolerance 1E-16.
        eps <- 1e-16
        if(X2mcd@crit < obj + eps)
        {
            X2cov <- getCov(X2mcd)
            X2center <- getCenter(X2mcd)
            if(trace)
                cat("\nFinal step - PC of MCD cov used.\n")
        }else
        {
            consistencyfactor <- median(mah)/qchisq(0.5,k)
            mah <- mah/consistencyfactor
            weights <- ifelse(mah <= qchisq(0.975, k), TRUE, FALSE)

##          VT::27.08.2010 - not necessary, cov.wt is doing it properly
##            wcov <- .wcov(X2, weights)
            wcov <- cov.wt(x=X2, wt=weights, method="ML")
            X2center <- wcov$center
            X2cov <- wcov$cov
            if(trace)
                cat("\nFinal step - PC of a reweighted cov used.\n")
        }

        ee <- eigen(X2cov)
        P6 <- ee$vectors

        center <- as.vector(center + X2center %*% t(rot))
        eigenvalues <- ee$values
        loadings <- rot %*% P6
##        scores <- (X2 - repmat(X2center, n, 1)) %*% P6
        scores <- (X2 - matrix(rep(X2center, times=n), nrow=n, byrow=TRUE)) %*% P6

        if(is.list(dimnames(data)) && !is.null(dimnames(data)[[1]]))
        {
            dimnames(scores)[[1]] <- dimnames(data)[[1]]
        } else {
            dimnames(scores)[[1]] <- 1:n
        }
        dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))
        dimnames(scores)[[2]] <- as.list(paste("PC", seq_len(ncol(scores)), sep = ""))

        res <- new("PcaHubert",call=cl,
                            loadings=loadings,
                            eigenvalues=eigenvalues,
                            center=center,
                            scale=myscale,
                            scores=scores,
                            k=k,
                            quan=h,
                            alpha=alpha,
                            n.obs=n)

    }

    cl[[1]] <- as.name("PcaHubert")
    res@call <- cl

    ## Compute distances and flags
    res <- pca.distances(res, data, Xsvd$rank, crit.pca.distances)
    return(res)
}


##
##  Returns 'ndirect' (random) directions through the data -
##      each through a pair of (randomly choosen) data points.
##  If all=TRUE all possible directions (all possible pairs of points)
##      are returned.
##
extradir <- function(data, ndirect, all=TRUE){
    if(all)                             # generate all possible directions (all pairs of n points)
    {
        cc <- combn(nrow(data), 2)
        B2 <- data[cc[1,],] - data[cc[2,],]
    }
    else {                              # generate 'ndirect' random directions
        if(TRUE){
                uniran <- function(seed = 0){
                    seed<-floor(seed*5761)+999
                    quot<-floor(seed/65536)
                    seed<-floor(seed)-floor(quot*65536)
                    random<-seed/65536
                    list(seed=seed, random=random)
                }

                ##
                ##  Draws a random subsubsample of k objects out of n.
                ##  This function is called if not all (p+1)-subsets
                ##  out of n will be considered.
                ##
                randomset <- function(n, k, seed){
                    ranset <- vector(mode="numeric", length=k)
                    for(j in 1:k){
                        r <- uniran(seed)
                        seed <- r$seed
                        num <- floor(r$random * n) + 1

                        if(j > 1){
                            while(any(ranset == num)){
                                r <- uniran(seed)
                                seed <- r$seed
                                num <- floor(r$random * n) + 1
                            }
                        }

                        ranset[j] <- num
                    }
                    ans<-list()
                    ans$seed <- seed
                    ans$ranset <- ranset
                    ans
                }

                n <- nrow(data)
                p <- ncol(data)
                r <- 1
                B2 <- matrix(0,ndirect, p)
                seed <- 0
                while(r <= ndirect) {
                    sseed <- randomset(n, 2, seed)
                    seed  <- sseed$seed
                    B2[r,] <- data[sseed$ranset[1], ] - data[sseed$ranset[2],]
                    r <- r + 1
                }
        } else
        {
                B2 <- matrix(0,ndirect, ncol(data))
                for(r in 1:ndirect) {
                    smpl <- sample(1:nrow(data), 2)                 # choose a random pair of points
                    B2[r,] <- data[smpl[1], ] - data[smpl[2], ]     # calculate the random direction based on these points
                }
        }

    }
    return(B2)
}

unimcd <- function(y, quan){
    out <- list()
    ncas <- length(y)
    len <- ncas-quan+1

    if(len == 1){
        out$tmcd <- mean(y)
        out$smcd <- sqrt(var(y))
    } else {
        ay <- c()
        I <- order(y)
        y <- y[I]
        ay[1] <- sum(y[1:quan])
        for(samp in 2:len){
            ay[samp]<-ay[samp-1]-y[samp-1]+y[samp+quan-1]
        }
        ay2<-ay^2/quan
        sq<-c()
        sq[1]<-sum(y[1:quan]^2)-ay2[1]
        for(samp in 2:len){
            sq[samp]<-sq[samp-1]-y[samp-1]^2+y[samp+quan-1]^2-ay2[samp]+ay2[samp-1]
        }
        sqmin<-min(sq)
        Isq<-order(sq)
        ndup<-sum(sq == sqmin)
        ii<-Isq[1:ndup]
        slutn<-c()
        slutn[1:ndup]<-ay[ii]
        initmean<-slutn[floor((ndup+1)/2)]/quan
        initcov<-sqmin/(quan-1)
        res<-(y-initmean)^2/initcov
        sortres<-sort(res)
        factor<-sortres[quan]/qchisq(quan/ncas,1)
        initcov<-factor*initcov
        res<-(y-initmean)^2/initcov
        quantile<-qchisq(0.975,1)
        out$weights<-(res<quantile)
        out$tmcd<-sum(y*out$weights)/sum(out$weights)
        out$smcd<-sqrt(sum((y-out$tmcd)^2*out$weights)/(sum(out$weights)-1))
        Iinv<-order(I)
        out$weights<-out$weights[Iinv]
    }
    return(out)
}
