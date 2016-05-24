ProcrustesRotation <- function(X,target,translate=FALSE,dilate=FALSE) {
        C <- crossprod(target,X)
        S <- svd(C)
        A <- S$v %*% t(S$u)
        if(dilate) A <- A * sum(S$d) / sum(diag(crossprod(t(X))))
        if(translate){
            centX <- apply(X,2,mean)
            centY <- apply(target,2,mean)
            b <- centY - crossprod(A,centX)
            f <- function(X) sweep(X %*% A,2,-b)
        } else f <- function(X) X %*% A
        f
}

procrustes <- function(x,...) UseMethod("procrustes")
procrustes.unfolding <- function(x,use=attr(x,"procrustes_use"),target,...){

    if(!length(use)) use <- "B"

    X <- switch(use, A = x$A, B = x$B)

    if(is.list(target)){

        nms <- unique(unlist(lapply(target,names)))
        if(length(nms)) nc <- length(nms)
        else nc <- max(unlist(lapply(target,length)))

        nr <- length(target)
        tmp <- matrix(0,nrow=nr,ncol=nc)

        if(length(lnms <- names(target)))
          rownames(tmp) <- lnms

        if(length(nms)){

          colnames(tmp) <- nms
          for(i in seq_along(target)){

            ti <- target[[i]]
            tmp[i,names(ti)] <- ti
          }
        }
        else {

          for(i in seq_along(target)){

            ti <- target[[i]]
            tmp[i,seq_along(ti)] <- ti
          }
        }
        target <- tmp
    }


    if(ncol(target)>ncol(X)) stop("undefined columns selected")

    if(length(colnames(target))) {

        rotdims <- match(colnames(target),colnames(X))
        if(any(is.na(rotdims))) stop("undefined columns selected")
        }
    else rotdims <- seq_len(ncol(target))

    tmp <- matrix(0,nrow=nrow(X),ncol=ncol(target))

    if(length(rownames(target))){

      ii <- match(rownames(target),rownames(X))
      if(any(is.na(ii))) stop("undefined rows selected")
    }
    else {

      ii <- seq_len(nrow(target))
    }

    tmp[ii,] <- target
    target <- tmp


    rot <- ProcrustesRotation(X[,rotdims,drop=FALSE],target,...)

    x$A[,rotdims] <- rot(x$A[,rotdims,drop=FALSE])
    x$B[,rotdims] <- rot(x$B[,rotdims,drop=FALSE])

    x
}