# shrinking clusters
shrinking <- function(y, K, disMethod = "Euclidean", 
  eps = 1.0e-4, itmax = 20)
{
    disMethod = match.arg(disMethod, c("Euclidean", "1-corr"))
 
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }
    y <- as.matrix(y)
    n <- nrow(y)
    p <- ncol(y)
 
    if(disMethod == "Euclidean") {
        disMethod2 <- 1
    } else { #if (disMethod == "corr") {
        disMethod2 <- 2
    } 
 
    K2 <- K + 1
 
    ynew <- matrix(0, nrow = n, ncol = p)
  
    res <- .Fortran("sharpen", 
        as.double(y), 
        as.integer(n), 
        as.integer(p), 
        as.integer(K2), 
        as.integer(K), 
        as.integer(itmax), 
        as.double(eps),
        as.integer(disMethod2), 
        ynew = as.double(ynew), 
        PACKAGE = "clues") 
 
    # Like R, Fortran stores a matrix by columns
    y.new <- matrix(res$ynew, nrow = n, ncol = p, byrow = FALSE)
    invisible(y.new)
}


