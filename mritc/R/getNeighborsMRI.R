getNeighborsMRI <- function(mask, nnei){
    if(length(dim(mask)) != 3)
        stop("The 'mask' has to be of dimension 3.")
    if(! all(unique(as.vector(mask)) %in% 0:1))
        stop("The value of mask has to be either 1 or 0.")
    if(all(mask==0))
        stop("All voxels are outside the mask.")
    if(! nnei %in% c(6, 18, 26))
        stop("'nnei' has to be among (6, 18, 26)")

    expand6 <- matrix(c(-1, 0, 0,
                        1, 0, 0,
                        0, -1, 0,
                        0, 1, 0,
                        0, 0, -1,
                        0, 0, 1), ncol=6)
    
    expand18 <- matrix(c(-1, 0, 0,
                         +1, 0, 0,
                         0,  -1, 0,
                         0,  +1, 0,
                         -1, -1, 0,
                         +1, +1, 0,
                         +1, -1, 0,
                         -1, +1, 0,
                         0, 0, -1,
                         0, 0, 1,
                         -1, 0, -1,
                         +1, 0, +1,
                         +1, 0, -1,
                         -1, 0, +1,
                         0, -1, -1,
                         0, +1, +1,
                         0, +1, -1,
                         0, -1, +1), ncol=18)
    expand26 <- matrix(c(-1, 0, 0,
                         +1, 0, 0,
                         0,  -1, 0,
                         0,  +1, 0,
                         -1, -1, 0,
                         +1, +1, 0,
                         +1, -1, 0,
                         -1, +1, 0,
                         0, 0, -1,
                         0, 0, 1,
                         -1, 0, -1,
                         +1, 0, +1,
                         +1, 0, -1,
                         -1, 0, +1,
                         0, -1, -1,
                         0, +1, +1,
                         0, +1, -1,
                         0, -1, +1,
                         -1, -1, -1,
                         -1, +1, +1,
                         -1, +1, -1,
                         -1, -1, +1,
                         +1, -1, -1,
                         +1, +1, +1,
                         +1, +1, -1,
                         +1, -1, +1), ncol=26)

    d <- dim(mask)
    n <- sum(mask)
    newmask <- array(0, dim=d+2)
    newmask[1:d[1]+1, 1:d[2]+1, 1:d[3]+1] <- mask
    mask <- newmask

    maskn <- array(cumsum(mask), dim=dim(mask))
    maskn[mask==0] <- n + 1

    focus <- which(mask==1, arr.ind=TRUE)
    focus <- t(focus)
    
    if(nnei == 6)
        neighbors <- do.call(cbind, lapply(1:6, function(i)
                         maskn[t(focus + expand6[,i])]))
    else{
        if(nnei == 18){
            neighbors <- do.call(cbind, lapply(1:18, function(i)
                         maskn[t(focus + expand18[,i])]))
        }
        else{
            neighbors <- do.call(cbind, lapply(1:26, function(i)
                         maskn[t(focus + expand26[,i])]))
        }
            
    }
    neighbors
}

