pre <- function(mask){
    n <- sum(mask==1)
    if(is.vector(mask)){
        mask <- c(0, mask, 0)
        mask.new <- cumsum(mask)
        mask.new[mask==0] <- n + 1
    }
    else if(is.matrix(mask)){
        nr <- nrow(mask)
        nc <- ncol(mask)
        mask.new0 <- matrix(0, nr + 2, nc + 2)
        mask.new0[2:(nr+1), 2:(nc+1)] <- mask
        mask.new <- cumsum(mask.new0)
        mask.new[mask.new0==0] <- n + 1
    }
    else{
        nr <- dim(mask)[1]
        nc <- dim(mask)[2]
        nz <- dim(mask)[3]
        mask.new0 <- array(0, dim=c(nr + 2, nc + 2, nz + 2))
        mask.new0[2:(nr+1), 2:(nc+1), 2:(nz + 1)] <- mask
        mask.new <- cumsum(mask.new0)
        mask.new[mask.new0==0] <- n + 1
    }
    list(mask=mask.new, n=n)
}

get2NeighborsLine <- function(mask, n){
    focus <- which(mask <= n)
    cbind(mask[focus-1], mask[focus+1])
}

get2NeighborsPlane <- function(mask, n, nr, neiStruc){
    neighbors <- vector("list", 4)
    focus <- which(mask <= n)
    if(neiStruc[1]==T)
        neighbors[[1]] <- cbind(mask[focus-1], mask[focus+1])
    if(neiStruc[2]==T)
        neighbors[[2]] <- cbind(mask[focus-nr], mask[focus+nr])
    if(neiStruc[3]==T)
        neighbors[[3]] <- cbind(mask[focus-nr-1], mask[focus+nr+1])
    if(neiStruc[4]==T)
        neighbors[[4]] <- cbind(mask[focus-nr+1], mask[focus+nr-1])
    neighbors
}

get2NeighborsCube <- function(mask, n, nr, nc, neiStruc){
    neighbors <- vector("list", 12)
    focus <- which(mask <= n)
    n12 <- nr*nc
    if(neiStruc[1]==T)
        neighbors[[1]] <- cbind(mask[focus-1], mask[focus+1])
    if(neiStruc[2]==T)
        neighbors[[2]] <- cbind(mask[focus-nr], mask[focus+nr])
    if(neiStruc[3]==T)
        neighbors[[3]] <- cbind(mask[focus-nr-1], mask[focus+nr+1])
    if(neiStruc[4]==T)
        neighbors[[4]] <- cbind(mask[focus-nr+1], mask[focus+nr-1])

    if(neiStruc[5]==T)
        neighbors[[5]] <- cbind(mask[focus-1], mask[focus+1])
    if(neiStruc[6]==T)
        neighbors[[6]] <- cbind(mask[focus-n12], mask[focus+n12])
    if(neiStruc[7]==T)
        neighbors[[7]] <- cbind(mask[focus-n12-1], mask[focus+n12+1])
    if(neiStruc[8]==T)
        neighbors[[8]] <- cbind(mask[focus-n12+1], mask[focus+n12-1])

    if(neiStruc[9]==T)
        neighbors[[9]] <- cbind(mask[focus-nr], mask[focus+nr])
    if(neiStruc[10]==T)
        neighbors[[10]] <- cbind(mask[focus-n12], mask[focus+n12])
    if(neiStruc[11]==T)
        neighbors[[11]] <- cbind(mask[focus-n12-nr], mask[focus+n12+nr])
    if(neiStruc[12]==T)
        neighbors[[12]] <- cbind(mask[focus-n12+nr], mask[focus+n12-nr])


    neighbors
}

get2NeighborsCube3DDiad <- function(mask){
    expand <- matrix(c(-1, -1, -1,
                       -1, +1, +1,
                       -1, +1, -1,
                       -1, -1, +1,
                       +1, -1, -1,
                       +1, +1, +1,
                       +1, +1, -1,
                       +1, -1, +1), ncol=8)
    browser()
    n <- sum(mask)
    maskn <- array(cumsum(mask), dim=dim(mask))
    maskn[mask==0] <- n + 1
    focus <- which(mask==1, arr.ind=TRUE)
    focus <- t(focus)
    neighbors <- do.call(cbind, lapply(1:8, function(i)
                         maskn[t(focus + expand[,i])])) 
    
}

getNeighborsB <- function(neiStruc, neighbors2){
    nvertex <- nrow(neighbors2)
    neighbors <- matrix(nvertex+1, nrow=nvertex+1, ncol=neiStruc)
    neighbors[1:nvertex,1:2] <- neighbors2
    if(neiStruc > 2)
        for(i in 1:(neiStruc/2-1)){
            neighbors[,i*2+1] <- neighbors[neighbors[,(i-1)*2+1],1]
            neighbors[,i*2+2] <- neighbors[neighbors[,(i-1)*2+2],2]
        }
    neighbors[-(nvertex+1),]
}

getNeighborsLine <- function(mask, neiStruc){
    r <- pre(mask)
    neighbors2 <- get2NeighborsLine(r$mask, r$n)
    getNeighborsB(neiStruc, neighbors2)

}

getNeighborsPlane <- function(mask, neiStruc){
    r <- pre(mask)
    nr <- nrow(mask) + 2
    neighbors2 <- get2NeighborsPlane(r$mask, r$n, nr, neiStruc > 0)
    if(all(neiStruc <= 2))
        neighbors <- do.call(cbind, neighbors2)
    else{
        neighbors <- NULL
        for(i in 1:4){
            if(!is.null(neighbors2[[i]])){
                new <- getNeighborsB(neiStruc[i], neighbors2[[i]])
                neighbors <- cbind(neighbors, new)
            }
        }
    }
    neighbors
}

getNeighborsCube <- function(mask, neiStruc){
    r <- pre(mask)
    nr <- dim(mask)[1] + 2
    nc <- dim(mask)[2] + 2
    neiStruc <- as.vector(t(neiStruc))
    neighbors2 <- get2NeighborsCube(r$mask, r$n, nr, nc, neiStruc > 0)
    if(all(neiStruc <= 2))
        neighbors <- do.call(cbind, neighbors2)
    else{
        neighbors <- NULL
        for(i in 1:12){
            if(!is.null(neighbors2[[i]])){
                new <- getNeighborsB(neiStruc[i], neighbors2[[i]])
                neighbors <- cbind(neighbors, new)
            }
        }
    }
    neighbors
}


getNeighbors <- function(mask, neiStruc){
    if(!(is.vector(mask) || is.matrix(mask) || length(dim(mask))==3))
        stop("The graph has to be in 1D, 2D, or 3D.")

    if(length(mask) < 3)
        stop("There are at least 3 vertices.")
    if(is.vector(mask))
        if(!is.vector(neiStruc) || length(neiStruc) != 1)
            stop("For a line, you need to tell the number of neighbors
                 on only one direction.")
    if(is.matrix(mask))
        if(!is.vector(neiStruc) || length(neiStruc)!=4)
            stop("For a plane, you need to tell the number of neighbors
                 on four directions.")
    if(length(dim(mask))==3)
        if(!is.matrix(neiStruc) || nrow(neiStruc)!=3 || ncol(neiStruc)!=4)
            stop("For a cube, you need to tell the number of neighbors
                 on four directions from all three perspective.")

    if(! all(neiStruc >= 0))
        stop("Number of neighbors have to be >= 0.")
    if(! all(neiStruc %% 2 == 0))
        stop("Number of neighbors have to be all even.")

    if(is.vector(mask))
        neighbors <- getNeighborsLine(mask, neiStruc)

    if(is.matrix(mask))
        neighbors <- getNeighborsPlane(mask, neiStruc)

    if(length(dim(mask))==3)
        neighbors <- getNeighborsCube(mask, neiStruc)

    neighbors
}
