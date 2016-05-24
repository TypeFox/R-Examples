# this is from the AnalyzeFMRI package
Sim.3D.GRF <- function(d, voxdim, sigma, ksize, mask = NULL, type = c("field", "max")) {

    ## simulates a GRF with covariance matrix sigma of dimension d with voxel dimensions voxdim
    ## if type = "max" then just the maximum of the field is returned
    ## if type = "filed" then just the filed AND the maximum of the field are returned

    if(!is.null(mask) && sum(d[1:3] == dim(mask)) < 3)  return("mask is wrong size")
    if(length(d) != 3 && length(d) != 4) return("array x should be 3D or 4D")
    if(length(d) == 3) {
        d <- c(d, 1)
        tmp <- 1
    }

    if((2 * floor(ksize / 2)) == ksize) stop(paste("ksize must be odd"))

    if(is.null(mask)) mask <- array(1, dim = d[1:3])
    space <- 1 + (type == "field") * prod(d)

    filtermat <- GaussSmoothKernel(voxdim, ksize, sigma)

    a <- .C("sim_grf",
            as.integer(d),
            as.double(aperm(filtermat, c(3, 2, 1))),
            as.integer(ksize),
            as.integer(aperm(mask, c(3, 2, 1))),
            as.integer((type == "field")),
            mat = double(space),
            max = double(1),
            PACKAGE = "neuRosim")

    if(type == "field") {
        mat <- array(a$mat, dim = d[4:1])
        mat <- aperm(mat, 4:1)
        if(tmp == 1) mat <- mat[, , , 1]
        return(list(mat = mat, max = a$max))
    }

    return(list(mat = NULL, max = a$max))


}

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


