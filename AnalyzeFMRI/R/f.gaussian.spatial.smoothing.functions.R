#functions to apply gaussian spatial smoothing to an array

    

GaussSmoothKernel<-function(voxdim = c(1 , 1, 1), ksize = 5, sigma = diag(3, 3))
#calculates a discretized smoothing kernel in up to 3 dimensions given an arbitrary covariance matrix
#sigma is covariance matrix of the gaussian
#doesn't have to be non-singular; zero on the diagonal of sigma indicate no smoothing in that direction
  
{
    if((2 * floor(ksize / 2)) == ksize) stop(paste("ksize must be odd"))
    
    a <- array(0, dim = c(ksize, ksize, ksize))
    centre <- (ksize + 1) / 2
    
    sig.ck <- c(TRUE, TRUE, TRUE)
    for(i in 1:3){
        if(sigma[i, i] == 0){
            sigma[i, i] <- 1
            sig.ck[i] <- FALSE
        }
    }
    sig.inv <- solve(sigma)
    sig.det <- abs(det(sigma))
    
    
    
    for(i in 1:ksize) {
        for(j in 1:ksize) {
            for(k in 1:ksize) {
                x <- (c(i, j, k) - centre) * voxdim
                a[i, j, k] <- ((2 * pi)^(-3 / 2)) * exp(-.5 * (t(x) %*% sig.inv %*% x)) / sqrt(sig.det)
            }
        }
    }
    if(sig.ck[1] == FALSE) a[-centre, , ] <- 0
    if(sig.ck[2] == FALSE) a[, -centre, ] <- 0
    if(sig.ck[3] == FALSE) a[, , -centre] <- 0
    a <- a / sum(a)

    return(a)
}


GaussSmoothArray <- function(x, voxdim = c(1, 1, 1), ksize = 5, sigma = diag(3, 3), mask = NULL, var.norm= FALSE )
{
    filtmat <- GaussSmoothKernel(voxdim, ksize, sigma)
    
    if(!is.array(x)) return("x should be an array")
    if(length(dim(x)) != 3 && length(dim(x)) != 4) return("array x should be 3D or 4D")
    tmp <- FALSE
    if(length(dim(x)) == 3) {
        x <- array(x, dim = c(dim(x), 1))
        tmp <- TRUE
    }
    if(is.null(mask)) mask <- array(1, dim = dim(x)[1:3])
    
    if(var.norm){
        d <- .Fortran("gaussfilter2",
                      as.double(x),
                      as.integer(dim(x)[1]),
                      as.integer(dim(x)[2]),
                      as.integer(dim(x)[3]),
                      as.integer(dim(x)[4]),
                      as.double(filtmat),
                      as.integer(ksize),
                      as.double(mask),
                      double(length(x)),
                      PACKAGE = "AnalyzeFMRI")
        c1 <- array(d[[9]], dim = dim(x))
        if(tmp) c1 <- c1[, , , 1]
    }
        
    else
    {   d <- .Fortran("gaussfilter1",
                      as.double(x),
                      as.integer(dim(x)[1]),
                      as.integer(dim(x)[2]),
                      as.integer(dim(x)[3]),
                      as.integer(dim(x)[4]),
                      as.double(filtmat),
                      as.integer(ksize),
                      as.double(mask),
                      as.double(x),
                      PACKAGE = "AnalyzeFMRI")
        c1 <- array(d[[9]], dim = dim(x))
        if(tmp) c1 <- c1[, , , 1]
    }
    
    return(c1)
}

    
