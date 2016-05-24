logLik.mpp <- function(object, SNOWcluster=NULL, ...){
    params <- object$params
    data <- object$data
    cif <- object$gif
    gparams <- eval(object$gmap)
    mparams <- eval(object$mmap)
    TT <- object$TT
    #---------------------------------------------
    #   evaluate  sum(log(lambda(t)))
    evalpts0 <- data[((data[,"time"] <= TT[2]) &
                     (data[,"time"] >= TT[1])), ]
    if (!is.null(SNOWcluster)){
        if (!requireNamespace("parallel")) stop("The R package 'parallel' is required to use the argument SNOWcluster")
        if (inherits(SNOWcluster, "cluster")){
            #  swap order of evalpts & data, TT & tplus not needed
            gif_cluster <- function(evalpts, data, params)
                cif(data, evalpts, params, TT=NA, tplus=FALSE)
            n <- nrow(evalpts0)
            m <- length(SNOWcluster)
            #  adjustment for different CPU speeds
            if (is.null(attr(SNOWcluster, "cpu.spd"))) cpu.spd <- rep(1, m)
            else cpu.spd <- attr(SNOWcluster, "cpu.spd")
            cpu.spd <- cpu.spd/sum(cpu.spd)
            N <- n*(n+1)/2
#           n1 <- round(sqrt(2*N/m))
            n1 <- round(sqrt(2*N*cpu.spd[1]))
            w <- list()
            w[[1]] <- evalpts0[1:n1,]
            if (m > 2){
                for (i in 2:(m-1)){
#                   n2 <- round(sqrt((n1+0.5)^2 - 2*(n1-N/m)))
                    n2 <- round(sqrt((n1+0.5)^2 -
                                2*(n1-N*cpu.spd[i])))
                    w[[i]] <- evalpts0[(n1+1):n2,]
                    n1 <- n2
                }
            }
            w[[m]] <- evalpts0[(n1+1):n,]
            L0 <- parallel::clusterApply(SNOWcluster, w, gif_cluster,
                               data=data, params=gparams)
            L1 <- sum(log(L0[[1]]))
            for (i in 2:m) L1 <- L1 + sum(log(L0[[i]]))
        } else stop("ERROR: Object SNOWcluster has wrong class")
    } else L1 <- sum(log(cif(data, evalpts0, gparams)))
    #---------------------------------------------
    L2 <- sum(cif(data, NULL, gparams, TT=TT))
    dmarks <- object$marks[[1]]
    if (is.null(dmarks))
        L3 <- 0
    else
        L3 <- sum(dmarks(data, data, mparams))
    LL <- L1 - L2 + L3
    return(LL)
}

