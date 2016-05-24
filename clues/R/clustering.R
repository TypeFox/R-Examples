# clustering
clustering <- function(y, disMethod = "Euclidean")
{
    disMethod <- match.arg(disMethod, c("Euclidean", "1-corr"))
 
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }
    n <- as.integer(nrow(y))
    p <- as.integer(ncol(y))
 
    n1 <- n - 1
 
    if(disMethod == "Euclidean") {
        disMethod2 <- 1
    } else { 
        disMethod2 <- 2
    } 
 
    #output
    point <- rep(0, n)
    db <- rep(0.0, n)
    omin <- 0.0
    nClusters <- 0
    mem <- rep(0.0, n)
    size <- rep(0.0, n)
    
    res <- .Fortran("clustering", 
        as.double(y), 
        as.integer(n), 
        as.integer(n1), 
        as.integer(p),
        as.integer(disMethod2),  
        point = as.integer(point), 
        db = as.double(db), 
        omin = as.double(omin), 
        nClusters = as.integer(nClusters), 
        mem = as.integer(mem), 
        size = as.integer(size), 
        PACKAGE = "clues") 
 
    g <- res$nClusters
    size <- res$size[1:g]
    db <- res$db[1:n1]
    resList <- list(mem = res$mem, size = size, g = g, 
        db = db, point = res$point, omin = res$omin)
 
    return(resList)
}

