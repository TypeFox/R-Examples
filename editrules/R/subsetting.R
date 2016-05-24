
#' Row index operator for \code{editmatrix}
#'
#' Use this operator to select edits from an editmatrix or editarray object.
#'
#' @method [ editmatrix 
#' @param x an object of class \code{\link{editmatrix}} or \code{\link{editarray}}
#' @param i the row index in the edit matrix (numeric, logical or rowname)
#' @param j the column index in the edit matrix
#' @param ... arguments to be passed to other methods. Currently ignored.
#' @rdname subsetting
#' @export
`[.editmatrix` <- function(x, i, j, ...){
    if (!missing(i) && is.character(i) ) i <- match(i,rownames(x),nomatch=0)
    E <- neweditmatrix(
        A = as.matrix(x)[i, j, drop=FALSE],
        ops = getOps(x)[i]
    )
    attr(E,"H") <- attr(x,"H")[i, , drop=FALSE]
    attr(E,"h") <- attr(x,"h")
    E
}


#' Row index operator for \code{cateditmatrix}
#'
#' @method [ editarray
#' @rdname subsetting
#' @keywords internal
#'
`[.cateditmatrix` <- function(x, i, j, ...){
    if (!missing(i) && is.character(i) ) i <- match(i,rownames(x),nomatch=0)
    E <- as.editmatrix( 
        getA(x)[i,,drop=FALSE], 
        getb(x)[i],
        getOps(x)[i],
        binvars=attr(x,'binvars')
    )
    class(E) <- c("cateditmatrix","editmatrix")
    E
}

#' Row index operator for \code{editarray}
#'
#' @method [ editarray
#' @rdname subsetting
#' @export
#'
`[.editarray` <- function(x, i, j, ...){
    A <- getArr(x)[i,j,drop=FALSE]
    sep <- getSep(x)
    ind <- indFromArray(A,sep)
    H <- getH(x)
    if (!is.null(H)) H <- H[i,j,drop=FALSE]
    neweditarray(E=A, ind=ind, sep=sep, names=getnames(x)[i],levels=getlevels(x)[j],H=H)
}

#' Row index operator for \code{editset}
#' Note: the datamodel is not changed
#' @method [ editset
#' @rdname subsetting
#' @export
`[.editset` <- function(x,i,j, ...){
    if (is.numeric(i) && i[1] < 0){
      N <- nrow(x$num) + nrow(x$mixcat)
      # make a negative selection positive
      i <- seq_len(N)[i]
      print(list(i=i, N=N))
    }
    if ( is.logical(i) ) i <- which(i)
    nnum <- nrow(x$num)
    mixcat <- x$mixcat[i[i>nnum]-nnum]
    # remove edits from mixnum not occuring in mixcat
    v <- getVars(reduce(mixcat))
    mixnum <- x$mixnum[rownames(x$mixnum) %in% v,]
    # remove dummy variables from mixcat not referring to numerical edits anymore
    v <- getVars(x,type='dummy')
    delvars <- v[ !v %in% rownames(mixnum)]
    if ( length(delvars) > 0 ){
        ind <- getInd(mixcat)
        delcols <- do.call('c',ind[delvars])
        Amixcat <- getArr(mixcat)[,-delcols,drop=FALSE]
        sep <- getSep(mixcat)
        ind <- indFromArray(Amixcat,sep=sep)
        mixcat <- neweditarray(Amixcat,ind,sep) 
    }
    neweditset(
        num = x$num[i[i<=nnum]],
        mixnum = mixnum,
        mixcat = mixcat,
        condition = attr(x,"condition")
    )
}


#' Index operator for \code{editlist}
#' @method [ editlist
#' @rdname subsetting
#' @export
`[.editlist` <- function(x,i,j, ...){
    x <- unclass(x)[i]
    class(x) <- 'editlist'
    x
}



