maximum.dist <- function(data.x, data.y=data.x, rank=FALSE){
    dx <- dim(data.x)
    dy <- dim(data.y)
    
    if(is.null(dx) & is.null(dy) ){
        x.lab <- names(data.x)
        y.lab <- names(data.y)
        if(rank){
            nx <- length(data.x)
            ny <- length(data.y)
            if((nx==ny) && all(data.x==data.y)){
                    rx <- rank(data.x, na.last="keep", ties.method="average")/(nx+1)
                    ry <- rx
            }
            else{
                rxy <- rank(c(data.x, data.y), na.last="keep", ties.method="average")/(nx + ny+1)
                rx <- rxy[1:nx]
                ry <- rxy[-(1:nx)]
            }
            mdist <- abs(outer(rx, ry, FUN="-"))
        }
        else mdist <- abs(outer(data.x, data.y, FUN="-"))
    }
    else{
        if(is.null(dx) & !is.null(dy)){
            data.y <- data.matrix(data.y)
            if(is.list(data.x)) data.x <- data.matrix(data.x)
            else  data.x <- matrix(data.x, nrow=1)
        }
        if(!is.null(dx) & is.null(dy)){
            data.x <- data.matrix(data.x)
            if(is.list(data.y)) data.y <- data.matrix(data.y)
            else  data.y <- matrix(data.y, nrow=1)
        }
        else {
            data.x <- data.matrix(data.x)
            data.y <- data.matrix(data.y)
        }
        x.lab <- rownames(data.x)
        y.lab <- rownames(data.y)
        nx <- nrow(data.x)
        ny <- nrow(data.y)
        if(rank){
            if(all(dx==dy)){
                if(all(data.x==data.y)){
                    rx <- apply(data.x, 2, rank, na.last="keep", ties.method="average")/(nx+1)
                    ry <- rx
                }
                else {
                    rxy <- apply(rbind(data.x,data.y), 2, rank, na.last="keep", ties.method="average")/(nx+ny+1)
                    rx <- rxy[1:nx,, drop=FALSE]
                    ry <- rxy[-(1:nx), , drop=FALSE]
                }
            }
            else{
                rxy <- apply(rbind(data.x,data.y), 2, rank, na.last="keep", ties.method="average")/(nx+ny+1)
                rx <- rxy[1:nx, , drop=FALSE]
                ry <- rxy[-(1:nx), , drop=FALSE]
            }
        }
        else{
            rx <- data.x
            ry <- data.y
        }
        mdist <- matrix(0, nx, ny)
        for(i in 1:nx){
            dd <- abs(rx[i,] - t(ry))
            mdist[i,] <- apply(dd, 2, max, na.rm=TRUE)
        }
    }
    dimnames(mdist) <- list(x.lab, y.lab)
    mdist
}
