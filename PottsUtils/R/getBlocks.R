getBlocksLine <- function(mask, nblock){
    nvertex <- sum(mask == 1)
    lapply(1:nblock, function(x) seq(x, nvertex, by=nblock))
}

getBlocksPlane <- function(mask, nblock){
    ind <- which(mask == 1, arr.ind=T)
    oe1 <- ind[,1] %% 2
    oe2 <- ind[,2] %% 2
    if(nblock==2){
        mark <- abs(oe1 - oe2)
        lapply(1:2, function(x) which(mark == (x-1)))
    }
    else
        lapply(1:4, function(x) which((oe1 + oe2*2) == (x-1)))
}

getBlocksCube <- function(mask, nblock){
    ind <- which(mask == 1, arr.ind=T)
    x.mark <- ind[,1] %% 2
    y.mark <- ind[,2] %% 2
    z.mark <- ind[,3] %% 2

    if(nblock==2){
        BW <- ifelse(rowSums(cbind(x.mark, y.mark, z.mark)) %% 2 == 1, 1, 0)
        lapply(1:2, function(x) which(BW == (x-1)))
    }

    else
        lapply(1:8, function(x) which((x.mark + y.mark*2 + z.mark*4) == (x-1)))
}

getBlocks <- function(mask, nblock){
    if(!(is.vector(mask) || is.matrix(mask) || length(dim(mask))==3))
        stop("The graph has to be in 1D, 2D, or 3D.")
    if(nblock <= 0 || floor(nblock) != nblock)
        stop("The number of blocks has to be and integer bigger than 0")
    if(length(mask==1) < nblock)
        stop("The number of blocks has to be less than
             the number of vertices inside the mask.")
    if(is.matrix(mask) && !nblock %in% c(2, 4))
        stop("For a plane, the number of blocks has to be equal to 2 or 4.")
    if(length(dim(mask))==3 && !nblock %in% c(2, 8))
        stop("For a cube, the number of blocks has to be equal to 2 or 8.")

    if(is.vector(mask))
        getBlocksLine(mask, nblock)

    else if(is.matrix(mask))
        getBlocksPlane(mask, nblock)

    else if(length(dim(mask))==3)
        getBlocksCube(mask, nblock)
}



