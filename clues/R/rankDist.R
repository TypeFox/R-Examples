rankCorr <- function(dat, method) 
{
    res <- cor(t(dat), method = method)
    return(res)
}

rankDist <- function(dat, method = "spearman")
{
    if(!is.matrix(dat))
    { dat <- matrix(dat, ncol = 1) }
    mat <- rankCorr(dat, method)
    dMat <- as.dist((1 - mat) / 2)
    return(dMat)
}

  
getDisti <- function(i, dat, method = "euclidean")
{
    if(!is.matrix(dat))
    { dat <- matrix(dat, ncol = 1) }
    n <- nrow(dat)
    p <- ncol(dat)
    # not like R, C stores a matrix by row. 
    # Hence need to transform 'dat' before pass it to C
    dat2 <- as.vector(t(dat))
 
    tt <- .C("getDisti", 
        as.integer(i), 
        as.double(dat2),
        as.integer(n), 
        as.integer(p), 
        di = as.double(rep(0, n)), 
        PACKAGE = "clues")
 
    di <- tt$di
    return(di)
}

getDistij <- function(i, j, dat, method = "euclidean")
{
    if(!is.matrix(dat))
    { dat <- matrix(dat, ncol = 1) }
    p <- ncol(dat)
 
    # not like R, C stores a matrix by row. 
    # Hence need to transform 'dat' before pass it to C
    dat2 <- as.vector(t(dat))
    tt <- .C("getDistij", 
        as.integer(i), 
        as.integer(j), 
        as.double(dat2), 
        as.integer(p), 
        dij = as.double(0), 
        PACKAGE = "clues")
 
    dij <- tt$dij
    return(dij)
}

getDistSets <- function(set1, set2, dat, method = "euclidean")
{
    if(!is.matrix(dat))
    { dat <- matrix(dat, ncol = 1) }
    n1 <- length(set1)
    n2 <- length(set2)
    p <- ncol(dat)
 
    # not like R, C stores a matrix by row. 
    # Hence need to transform 'dat' before pass it to C
    dat2 <- as.vector(t(dat))
    tt <- .C("getDistSets", 
        as.integer(set1), 
        as.integer(n1), 
        as.integer(set2), 
        as.integer(n2), 
        as.double(dat2), 
        as.integer(p), 
        dMat = as.double(rep(0, n1 * n2)), 
        PACKAGE = "clues")
 
    # not like R, C stores a matrix by row. 
    # Hence need to set 'byrow=T'
    dMat <- matrix(tt$dMat, nrow = n1, ncol = n2, byrow = T)
    return(dMat)
}

