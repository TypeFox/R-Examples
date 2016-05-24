##  Robust SEQIMPUTE is a robust sequential imputation method.
##
##  The seqimpute method is described in:
##     S. Verboven, K. Vanden Branden and P. Goos (2007)
##     "Sequential Imputation for missing values", Computational Biology and
##     Chemistry, 31: 320-327.
##
##  The robust seqimpute method is described in:
##
##      .......
##
##  a) If the data set is complete, just return it without doing anything.
##  b) The algorithm needs at least p complete cases - if there are no or less than p
##      complete cases, find a way to impute them in some other way.
##
##  library(rrcov)
##  data(phosphor)
##  x <- phosphor[,1:2]
##  x[10,2] <- NA
##  x[15,1] <- NA
##  y <- impSeqRob(x)
##  y[10,2]     # 37.51664
##  y[15,1]     # 26.38606
##  cbind(phosphor[,1:2],x,y$x)
##

impSeqRob <- function(x, alpha=0.9){

    if(is.data.frame(x))
    {
        x <- data.matrix(x)
    }else if(!is.matrix(x))
    {
        x <- matrix(x, length(x), 1, dimnames = list(names(x), deparse(substitute(x))))
    }
    xcall <- match.call()

    n <- nrow(x)
    p <- ncol(x)
    isnanx = is.na(x) + 0
    risnanx = apply(isnanx,1,sum) #observations with missing values

    ntotmiss <- length(which(risnanx > 0))  # all missing
    ntotcompl <- n - ntotmiss               # all complete

    if(ntotmiss == 0)                   # if no missing data - calculate
    {
        outx <- .outlSD(x, alpha=alpha) # outlyingness of the complete
        return(outx)                    # matrix and return x

    }

    ##  sort according to percentage of missing values (so that first the observations
    ##  with smallest number of missing values are handled)

    ##[sortx,Ix] = sort(risnanx);
    sortx <- sort.int(risnanx, index.return=TRUE)
    sorth <- sort.int(sortx$ix, index.return=TRUE)
    x = x[sortx$ix,]

##isnanx = isnanx(Ix,:);
##risnanx = sortx; %observations with missings
    isnanx = is.na(x) + 0

    risnanx = sortx$x #observations with missing values

    complobs <- which(risnanx == 0)
    ncomplobs <- length(complobs)
    initncomplobs <- ncomplobs

    h <- h.alpha.n(alpha, ncomplobs, p)     # h <- floor(alpha*ncomplobs)
    misobs <- which(risnanx != 0)
    nmisobs <- length(misobs)

    mincomplete <- ceiling(p/alpha)         # minimum required complete observations
    if(ntotcompl >= 5*p)                    # more than enough complete - do nothing
    {
        ## do nothing
        r <- p
    }else if(ntotcompl >= mincomplete)      # if enough complete observations check
    {                                       # if the cov matrix is not singular
        covm <- .covSD(x[complobs,], h)
        r <- rankMM(covm$cov)
    }

    if(ntotcompl < mincomplete || r < p)        # if not enough complete observations
    {
        ## impute the missing data using package norm
        s <- prelim.norm(x)                     # do preliminary manipulations
        thetahat <- em.norm(s, showits=FALSE)   # find the mle
        rngseed(1234567)                        # set random number generator seed
        ximp <- imp.norm(s, thetahat, x)        # impute missing data under the MLE
        xx<-imp.norm(s, thetahat, x)            # impute missing data under the MLE

        outx <- .outlSD(ximp, alpha=alpha)      # outlyingness of the complete
        return(outx)                            # matrix and return x

        repeat
        {
            x[(ncomplobs+1):mincomplete,] <- ximp[(ncomplobs+1):mincomplete,]
            risnanx[(ncomplobs+1):mincomplete] <- 0

            complobs <- which(risnanx == 0)
            ncomplobs <- length(complobs)
            ## h <- floor(alpha*ncomplobs)
            h <- h.alpha.n(alpha, ncomplobs, p)
            covm <- .covSD(x[complobs,], h)
            r <- rankMM(covm$cov)
            if(r >= p)
                break

            cat("\nMin: ",mincomplete, " Rank: ",r, "Complete: ",ncomplobs, "\n")
            mincomplete <- mincomplete + 1
        }

        ncomplobs <- length(complobs)
        initncomplobs <- ncomplobs
        ## h <- floor(alpha*ncomplobs)
        h <- h.alpha.n(alpha, ncomplobs, p)
        misobs <- which(risnanx != 0)
        nmisobs <- length(misobs)
    }

    ##  use outlyingness to obtain an initial estimate for cov and mean of the
    ##  initial complete data set

##    nrich1 <- ncomplobs*(ncomplobs-1)/2
##    ndirect <- min(500,nrich1)  # note 500 instead of 250
##    true <- ndirect == nrich1
##    B <- .extradir(x[complobs,], ndirect, true)    # n*ri

##    Bnorm <- vector(mode="numeric", length=nrow(B))   # ndir x 1
##    Bnorm <- apply(B, 1, vecnorm)   #
##    Bnormr <- Bnorm[Bnorm > 1.E-12]         # choose only the nonzero length vectors
##    m <- length(Bnormr)                     # the number of directions that will be used
##    B <- as.matrix(B[Bnorm > 1.E-12,])      # B[m x ncol(X)]
##    A <- diag(1/Bnormr) %*% B               # A[m x ncol(X)]


    ## projected points in columns

##    Y <- as.matrix(x[complobs,]) %*% t(A)     # n x m - projections of the n points on each of the m directions
##    Z <- matrix(0,ncomplobs, m)               # n x m - to collect the outlyingness of each point on each direction

##    medY <- apply(Y, 2, median)
##    madY <- apply(Y, 2, mad)
##    Z <- abs(Y - matrix(1,ncomplobs,1) %*% medY) / (matrix(1,ncomplobs,1) %*% madY)

##    d <- apply(Z,1,max)
##    ds <- sort.int(d, index.return=TRUE)

    ## initial estimates for mean and cov
##    covx <- cov(x[complobs[ds$ix[1:h]],])
##    mx <- colMeans(x[complobs[ds$ix[1:h]],])

    covm <- .covSD(x[complobs,], h)
    covx <- covm$cov
    mx <- covm$mean
    A <- covm$A
    medY <- covm$medY
    madY <- covm$madY
    d <- covm$d
    ds <- sort.int(d, index.return=TRUE)

    flag <- matrix(1, 1, n)
    flag[complobs[ds$ix[1:h]]] <- 0
    outl <- matrix(0,n,1)
    outl[complobs] = d

    ncomplobsupdate  <- h
    ## start sequential imputation
    for(inn in 1:nmisobs)
    {
        mvar = as.logical(isnanx[misobs[inn],])
        xo = x[misobs[inn],!mvar]

        ## FIRST estimate missing part of x
        if(inn == 1)
        {
            icovx <- solve(covx)
        }else if(flag[misobs[inn-1]] == 0) # update inv(covx) if new point has entered
        {
            f <- 1/sqrt(ncomplobsupdate-1) * (x[misobs[inn-1],] - mxo)
            icovx <- .update_invcov((ncomplobsupdate-1)/(ncomplobsupdate-2)*icovx, f)
        }
        x[misobs[inn], mvar] <-
            mx[mvar] - solve(icovx[mvar,mvar]) %*% icovx[mvar, !mvar] %*% as.matrix(xo - mx[!mvar])

        ## THEN calculate the outlyingness
        Y <- x[misobs[inn], ] %*% t(A)  # 1*ndirect
        Zx <- abs(Y - medY) / madY

        ## ??? dx <- max(Zx, [], 2)';
        dx <- apply(Zx, 1 ,max)

        if(dx <= ds$x[h])  # not an outlier
        {
            flag[misobs[inn]] <- 0
            ncomplobsupdate <- ncomplobsupdate + 1
            h <- h+1
            mxo <- mx

            ## UPDATE MEAN AND COV
            mx <- ((ncomplobsupdate-1) * mxo + x[misobs[inn], ]) / ncomplobsupdate
            covx <- (ncomplobsupdate-2)/(ncomplobsupdate-1)*covx +
                1/(ncomplobsupdate-1) * tcrossprod(x[misobs[inn], ] - mxo, x[misobs[inn], ] - mxo)
        }

        outl[misobs[inn]] <- dx

        ## ds <- sort([ds,dx]);
        d <- c(d, dx)
        ds <- sort.int(d, index.return=TRUE)
        ncomplobs <- ncomplobs + 1
    }

    ## put the observations in the original order
    xseq <- x[sorth$ix, ]
    outl <- outl[sorth$ix]
    flag <- flag[sorth$ix]
    list(x=xseq, outl=outl, flag=flag)
}

##
##  Returns 'ndirect' (random) directions through the data -
##      each through a pair of (randomly choosen) data points.
##  If all=TRUE all possible directions (all possible pairs of points)
##      are returned.
##
.extradir <- function(data, ndirect, all=TRUE){
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

.update_invcov <- function(ainv, f)
{
    f <- as.vector(f)
    ff <- f %*% t(f)
    ainv - ainv %*% ff %*% ainv / as.numeric((1 + f %*% ainv %*% f))
}

###
##  Calculate the outlyingness of the observations in x (complete)
##      Used when the data matrix contains no missing data
##
.outlSD <- function(x, alpha=0.9)
{

    if(is.data.frame(x))
    {
        x <- data.matrix(x)
    }else if(!is.matrix(x))
    {
        x <- matrix(x, length(x), 1, dimnames = list(names(x), deparse(substitute(x))))
    }

    n <- nrow(x)
    p <- ncol(x)
    h <- floor(alpha*n)

    nrich1 <- n * (n-1)/2
    ndirect <- min(500,nrich1)  # note 500 instead of 250
    true <- ndirect == nrich1

    B <- .extradir(x, ndirect, true)    # n*ri

    Bnorm <- vector(mode="numeric", length=nrow(B))   # ndir x 1
    Bnorm <- apply(B, 1, vecnorm)   #
    Bnormr <- Bnorm[Bnorm > 1.E-12]         # choose only the nonzero length vectors
    m <- length(Bnormr)                     # the number of directions that will be used
    B <- as.matrix(B[Bnorm > 1.E-12,])      # B[m x ncol(X)]
    A <- diag(1/Bnormr) %*% B               # A[m x ncol(X)]


    ## projected points in columns

    Y <- x %*% t(A)                 # n x m - projections of the n points on each of the m directions
    Z <- matrix(0,n, m)        # n x m - to collect the outlyingness of each point on each direction

    medY <- apply(Y, 2, median)
    madY <- apply(Y, 2, mad)
    Z <- abs(Y - matrix(1,n,1) %*% medY) / (matrix(1,n,1) %*% madY)

    d <- apply(Z,1,max)
    ds <- sort.int(d, index.return=TRUE)

    ## initial estimates for mean and cov
    covx <- cov(x[ds$ix[1:h], ])
    mx <- colMeans(x[ds$ix[1:h], ])

    flag <- matrix(1,1,n)
    flag[ds$ix[1:h]] <- 0
    outl <- matrix(0,n,1)
    outl = d

    list(xseq = x, outl=as.vector(outl), flag=as.vector(flag))
}

##  use outlyingness to obtain an initial estimate for cov and mean of the
##  initial complete data set
.covSD <- function(x, h)
{
    n <- nrow(x)
    p <- ncol(x)

    nrich1 <- n*(n-1)/2
    ndirect <- min(500,nrich1)      # note 500 instead of 250
    true <- ndirect == nrich1
    B <- .extradir(x, ndirect, true)    # n*ri

    Bnorm <- vector(mode="numeric", length=nrow(B))   # ndir x 1
    Bnorm <- apply(B, 1, vecnorm)   #
    Bnormr <- Bnorm[Bnorm > 1.E-12]         # choose only the nonzero length vectors
    m <- length(Bnormr)                     # the number of directions that will be used
    B <- as.matrix(B[Bnorm > 1.E-12,])      # B[m x ncol(X)]
    A <- diag(1/Bnormr) %*% B               # A[m x ncol(X)]


    ## projected points in columns
    Y <- as.matrix(x) %*% t(A)              # n x m - projections of the n points on each of the m directions
    Z <- matrix(0,n, m)                # n x m - to collect the outlyingness of each point on each direction

    medY <- apply(Y, 2, median)
    madY <- apply(Y, 2, mad)
    Z <- abs(Y - matrix(1,n,1) %*% medY) / (matrix(1,n,1) %*% madY)

    d <- apply(Z,1,max)
    ds <- sort.int(d, index.return=TRUE)

    ## initial estimates for mean and cov
    covx <- cov(x[ds$ix[1:h],])
    mx <- colMeans(x[ds$ix[1:h],])

    list(cov=covx, mean=mx, d=d, A=A, medY=medY, madY = madY)
}
