# using CH index as the measure of the strength of the clusters
clues_CH <- function(y, n0 = 5, alpha = 0.05, eps = 1.0e-4, itmax = 20,  
  K2.vec = NULL, CH2 = -3, disMethod = "Euclidean", quiet = FALSE)
{
    if(is.null(K2.vec))
    {
      K2.vec=n0
    }

    disMethod <- match.arg(disMethod, c("Euclidean", "1-corr"))
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }
    # first pass
    second <- FALSE
    res <- ChooseK_CH(y=y, y2 = y, n0=n0, alpha=alpha, 
      eps=eps, itmax=itmax, second=second, 
      K2.vec=K2.vec, CH2=CH2, 
        disMethod=disMethod, quiet=quiet)
      
    # check if we need second pass
    # if length(g.vec) = 1 or g.vec[pos-1]=g.vec[pos]+1=g.vec[pos+1]-1,
    # then do NOT need second pass
    len <- length(res$g.vec)
    if(len > 1) # may need second pass
    {   pos <- which(res$g.vec == res$g)
        len2 <- length(pos)
        len3 <- length(res$g.vec)
        pos1 <- pos[1]
        pos2 <- pos[len2]
        g <- res$g
        CH2 <- res$CH
        K12 <- res$K.vec[pos1]
        K21 <- res$K.vec[pos2]
        if(pos1 > 1)
        {   g1 <- res$g.vec[pos1 - 1]
            K11 <- res$K.vec[pos1 - 1]
        }
        if(pos2 < len3)
        {   K22 <- res$K.vec[pos2 + 1]
            g2 <- res$g.vec[pos2 + 1]
        }
        if(pos1 > 1 && g1 - g > 1) # need second pass
        {   K2.vec <- c(K11, K12)
            # second pass
            second <- TRUE
            res2 <- ChooseK_CH(y=y, y2 = res$y.old1, n0=n0, alpha=alpha, 
              eps=eps, itmax=itmax, second=second, K2.vec=K2.vec, 
              CH2=CH2, disMethod=disMethod, quiet=quiet)
            if(res2$myupdate) # update the final partition
            {   CH2 <- res2$CH 
                res <- res2
            }
        }
        if(pos2 < len3 && (g - g2) > 1) # need second pass
        {   K2.vec <- c(K21, K22)
            # second pass
            second <- TRUE
            res3 <- ChooseK_CH(y=y, y2 = res$y.old2, n0=n0, alpha=alpha, 
              eps=eps, itmax=itmax, second=second, K2.vec=K2.vec, 
              CH2=CH2, disMethod=disMethod, quiet=quiet)

            if(res3$myupdate) { res <- res3 }
        }
    }
    return(res)
}

# first pass
# if second = 1, then it is the first pass
# if second = 2, then it is the second pass
# y --- original data
# y2 --- data used to do shrinking and clustering
ChooseK_CH <- function(y, y2, n0 = 5, alpha = 0.05, eps = 1.0e-4, itmax = 20, 
    second = FALSE, K2.vec = NULL, CH2 = -3, disMethod = "Euclidean", 
    quiet = FALSE)
{
    if(is.null(K2.vec))
    {
      K2.vec=n0
    }

    disMethod <- match.arg(disMethod, c("Euclidean", "1-corr"))
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }
    # step 1
    y <- as.matrix(y)
    y2 <- as.matrix(y2)
    dat <- y
    dat2 <- y2
    nObs <- nrow(dat)
    nVars <- ncol(dat)
    nClusters0 <- n0
 
    if(disMethod == "Euclidean") {
        disMethod2 <- 1
    } else { 
        disMethod2 <- 2
    } 
 
    # input
    nObs1 <- nObs - 1
    ITMAX <- itmax
    nNeiVec2 <- K2.vec
 
    # output
    avgsFinal <- 0
    sFinal <- rep(0, nObs)
    memFinal <- rep(0, nObs)
    nClustersFinal <- 0
    clustSizeFinal <- rep(0, nObs)
    nNeiFinal <- 0
    nClustVec <- rep(0, nObs)
    nNeiVec <- rep(0, nObs)
    myt <- 0
    datold1 <- matrix(0, nrow = nrow(dat), ncol = ncol(dat))
    datold2 <- matrix(0, nrow = nrow(dat), ncol = ncol(dat))
    myupdate <- FALSE
    indexFlag <- FALSE
    
    res <- .Fortran("chooseK", 
        as.double(dat), 
        as.double(dat2), 
        as.integer(nObs), 
        as.integer(nObs1),
        as.integer(nVars), 
        as.integer(nClusters0), 
        as.double(alpha), 
        as.double(eps),
        as.integer(ITMAX), 
        as.logical(second), 
        as.integer(nNeiVec2), 
        as.double(CH2),  
        as.integer(disMethod2), 
        as.logical(indexFlag),
        as.logical(quiet),
        indFinal = as.double(avgsFinal), 
        sFinal = as.double(sFinal),
        memFinal = as.integer(memFinal), 
        nClustersFinal = as.integer(nClustersFinal), 
        clustSizeFinal = as.integer(clustSizeFinal), 
        nNeiFinal = as.integer(nNeiFinal), 
        nClustVec = as.integer(nClustVec), 
        nNeiVec = as.integer(nNeiVec), 
        myt = as.integer(myt), 
        datold1 = as.double(datold1), 
        datold2 = as.double(datold2), 
        myupdate = as.logical(myupdate),
        PACKAGE = "clues")  
       
       resList <- list(K = res$nNeiFinal, 
           size = res$clustSizeFinal[1:res$nClustersFinal],
           mem = res$memFinal, 
           g = res$nClustersFinal, 
           CH = res$indFinal,
           K.vec = res$nNeiVec[1:res$myt], 
           g.vec = res$nClustVec[1:res$myt],
           myupdate = res$myupdate, 
           y.old1 = matrix(res$datold1, nrow = nrow(dat), 
               ncol = ncol(dat), byrow = FALSE),
           y.old2 = matrix(res$datold2, nrow = nrow(dat), 
               ncol = ncol(dat), byrow = FALSE))
 
    return(resList)
}
 
get_CH <- function(y, mem, disMethod = "Euclidean")
{
    disMethod <- match.arg(disMethod, c("Euclidean", "1-corr"))
 
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }
    n <- nrow(y)
    p <- ncol(y)
 
    g <- length(unique(mem))
    size <- tapply(mem, mem, length)
 
    if(disMethod == "Euclidean") {
        disMethod2 <- 1
    } else { 
        disMethod2 <- 2
    }   
 
    res <- .Fortran("chindex", 
        as.double(y), 
        as.integer(n), 
        as.integer(p), 
        as.integer(mem), 
        as.integer(g), 
        as.integer(size), 
        as.integer(disMethod2), 
        CH = as.double(0),
        PACKAGE = "clues") 
 
    CH <- res$CH
 
    return(CH)
}


