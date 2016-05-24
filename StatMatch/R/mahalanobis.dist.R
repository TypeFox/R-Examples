mahalanobis.dist <- function(data.x, data.y=NULL, vc=NULL){

    xx <- as.matrix(data.x)
    if(is.null(data.y)) yy <- as.matrix(data.x)
    else yy <- as.matrix(data.y)

    if(is.null(vc)){
        if(is.null(data.y)) vc <- var(xx)
        else vc <- var(rbind(xx,yy))
    }

    ny <- nrow(yy)
    md <- matrix(0,nrow(xx), ny)
    for(i in 1:ny){
        md[,i] <- mahalanobis(xx, yy[i,], cov=vc)
    }
    if(is.null(data.y)) dimnames(md) <- list(rownames(data.x), rownames(data.x))
    else dimnames(md) <- list(rownames(data.x), rownames(data.y))
    sqrt(md)
}
