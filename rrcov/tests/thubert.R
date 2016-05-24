dodata <- function(nrep=1, time=FALSE, short=FALSE, full=TRUE, method=c("hubert", "hubert.mcd", "locantore", "cov", "classic")){
## Test the function PcaHubert() on the literature datasets:
##
## Call PcaHubert() for all regression datasets available in
##  robustbase/rrcov and print:
##  - execution time (if time == TRUE)
##  - loadings
##  - eigenvalues
##  - scores
##

    dopca <- function(x, xname, nrep=1){

        n <- dim(x)[1]
        p <- dim(x)[2]
        if(method == "hubert.mcd")
            pca <- PcaHubert(x, k=p)
        else if(method == "hubert")
            pca <- PcaHubert(x, mcd=FALSE)
        else if(method == "locantore")
            pca <- PcaLocantore(x)
        else if(method == "cov")
            pca <- PcaCov(x)
        else if(method == "classic")
            pca <- PcaClassic(x)
        else
            stop("Undefined PCA method: ", method)


        e1 <- getEigenvalues(pca)[1]
        e2 <- getEigenvalues(pca)[2]
        k <- pca@k

        if(time){
           xtime <- system.time(dorep(x, nrep, method))[1]/nrep
           xres <- sprintf("%3d %3d %3d %12.6f %12.6f %10.3f\n", dim(x)[1], dim(x)[2], k, e1, e2, xtime)
        }
        else{
            xres <- sprintf("%3d %3d %3d %12.6f %12.6f\n", dim(x)[1], dim(x)[2], k, e1, e2)
        }
        lpad<-lname-nchar(xname)
        cat(pad.right(xname, lpad), xres)

        if(!short){
            cat("Scores: \n")
            print(getScores(pca))

            if(full){
                cat("-------------\n")
                show(pca)
            }
            cat("----------------------------------------------------------\n")
        }
    }

    stopifnot(length(nrep) == 1, nrep >= 1)
    method <- match.arg(method)

    options(digits = 5)
    set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed

    lname <- 20

    ## VT::15.09.2013 - this will render the output independent
    ##  from the version of the package
    suppressPackageStartupMessages(library(rrcov))

    data(Animals, package = "MASS")
    brain <- Animals[c(1:24, 26:25, 27:28),]

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")

    cat("Data Set               n   p   k           e1           e2\n")
    cat("==========================================================\n")
    dopca(heart[, 1:2], data(heart), nrep)
    dopca(starsCYG, data(starsCYG), nrep)
    dopca(data.matrix(subset(phosphor, select = -plant)), data(phosphor), nrep)
    dopca(stack.x, data(stackloss), nrep)
    dopca(data.matrix(subset(coleman, select = -Y)), data(coleman), nrep)
    dopca(data.matrix(subset(salinity, select = -Y)), data(salinity), nrep)
##    dopca(data.matrix(subset(wood, select = -y)), data(wood), nrep)       # differences between the architectures
    dopca(data.matrix(subset(hbk,  select = -Y)),data(hbk), nrep)

##    dopca(brain, "Animals", nrep)
    dopca(milk, data(milk), nrep)
    dopca(bushfire, data(bushfire), nrep)
    cat("==========================================================\n")
}

dogen <- function(nrep=1, eps=0.49, method=c("hubert", "hubert.mcd", "locantore", "cov")){

    dopca <- function(x, nrep=1){
        gc()
        xtime <- system.time(dorep(x, nrep, method))[1]/nrep
        cat(sprintf("%6d %3d %10.2f\n", dim(x)[1], dim(x)[2], xtime))
        xtime
    }

    set.seed(1234)

    ## VT::15.09.2013 - this will render the output independent
    ##  from the version of the package
    suppressPackageStartupMessages(library(rrcov))
    library(MASS)

    method <- match.arg(method)

    ap <- c(2, 5, 10, 20, 30)
    an <- c(100, 500, 1000, 10000, 50000)

    tottime <- 0
    cat("     n   p       Time\n")
    cat("=====================\n")
    for(i in 1:length(an)) {
        for(j in 1:length(ap)) {
            n <- an[i]
            p <- ap[j]
            if(5*p <= n){
                xx <- gendata(n, p, eps)
                X <- xx$X
                ## print(dimnames(X))
                tottime <- tottime + dopca(X, nrep)
            }
        }
    }

    cat("=====================\n")
    cat("Total time: ", tottime*nrep, "\n")
}

dorep <- function(x, nrep=1, method=c("hubert", "hubert.mcd", "locantore", "cov")){

    method <- match.arg(method)
    for(i in 1:nrep)
    if(method == "hubert.mcd")
        PcaHubert(x)
    else if(method == "hubert")
        PcaHubert(x, mcd=FALSE)
    else if(method == "locantore")
        PcaLocantore(x)
    else if(method == "cov")
        PcaCov(x)
    else
        stop("Undefined PCA method: ", method)
}

#### gendata() ####
# Generates a location contaminated multivariate
# normal sample of n observations in p dimensions
#    (1-eps)*Np(0,Ip) + eps*Np(m,Ip)
# where
#    m = (b,b,...,b)
# Defaults: eps=0 and b=10
#
gendata <- function(n,p,eps=0,b=10){

    if(missing(n) || missing(p))
        stop("Please specify (n,p)")
    if(eps < 0 || eps >= 0.5)
        stop(message="eps must be in [0,0.5)")
    X <- mvrnorm(n,rep(0,p),diag(1,nrow=p,ncol=p))
    nbad <- as.integer(eps * n)
    xind <- vector("numeric")
    if(nbad > 0){
        Xbad <- mvrnorm(nbad,rep(b,p),diag(1,nrow=p,ncol=p))
        xind <- sample(n,nbad)
        X[xind,] <- Xbad
    }
    list(X=X, xind=xind)
}

pad.right <- function(z, pads)
{
### Pads spaces to right of text
    padding <- paste(rep(" ", pads), collapse = "")
    paste(z, padding, sep = "")
}

whatis<-function(x){
    if(is.data.frame(x))
        cat("Type: data.frame\n")
    else if(is.matrix(x))
        cat("Type: matrix\n")
    else if(is.vector(x))
        cat("Type: vector\n")
    else
        cat("Type: don't know\n")
}

#################################################################
##  VT::27.08.2010
##  bug report from Stephen Milborrow
##
test.case.1 <- function()
{
    X <- matrix(c(
          -0.79984, -1.00103,  0.899794,  0.00000,
           0.34279,  0.52832, -1.303783, -1.17670,
          -0.79984, -1.00103,  0.899794,  0.00000,
           0.34279,  0.52832, -1.303783, -1.17670,
           0.34279,  0.52832, -1.303783, -1.17670,
           1.48542,  0.66735,  0.716162,  1.17670,
          -0.79984, -1.00103,  0.899794,  0.00000,
           1.69317,  1.91864, -0.018363,  1.76505,
          -1.00759, -0.16684, -0.385626,  0.58835,
          -0.79984, -1.00103,  0.899794,  0.00000), ncol=4, byrow=TRUE)

    cc1 <- PcaHubert(X, k=3)

    cc2 <- PcaLocantore(X, k=3)
    cc3 <- PcaCov(X, k=3, cov.control=CovControlSest())

    cc4 <- PcaProj(X, k=2)           # with k=3 will produce warnings in .distances - too small eignevalues
    cc5 <- PcaGrid(X, k=2)           # dito

    list(cc1, cc2, cc3, cc4, cc5)
}

## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcov))

dodata(method="classic")
dodata(method="hubert.mcd")
dodata(method="hubert")

##dodata(method="locantore")
##dodata(method="cov")
test.case.1()
