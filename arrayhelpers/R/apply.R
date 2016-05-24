##' Conservative version of apply: tries to preserve the dimensions of the original array
##'
##' \code{\link[base]{apply}} returns its results as an array of shape (not-margin dimensions *
##' length of FUN's result) x MARGIN dimensions.
##'
##' \code{applycons} rearranges these results, so that
##'
##' 1. the result dimensions are appended rather than prepended, and
##' 2. tries to preserve the shape of the 
##' 
##' @param X an array or matrix.
##' @param MARGIN a vector giving the subscripts which the function will be applied over. E.g., for a
##' matrix \code{1} indicates rows, \code{2} indicates columns, \code{c(1, 2)} indicates rows and
##' columns. Where \code{X} has named dimnames, it can be a character vector selecting dimension
##' names.
##' @param FUN the function to be applied: see \sQuote{Details}. In the case of functions like
##' \code{+}, \code{\%*\%}, etc., the  function name must be backquoted or quoted.
##' @param ... optional arguments to \code{FUN}.
##' @return array
##' @author Claudia Beleites
##' @seealso \code{\link[base]{apply}}
##' @include arrayhelpers.R
##' @noRd
.applycons <- function (X, MARGIN, FUN, ...){
    FUN <- match.fun(FUN)

    ## Ensure that X is an array object
    dl <- length(dim(X))
    
    if(!dl) stop("dim(X) must have a positive length")

    if(is.object(X))
      X <- if(dl == 2L) as.matrix(X) else as.array(X)
    ## now record dim as coercion can change it
    ## (e.g. when a data frame contains a matrix).
    d <- dim(X)
    dn <- dimnames(X)
    ds <- seq_len(dl)

    ## Extract the margins and associated dimnames
    if (is.character(MARGIN)) {
        if(is.null(dnn <- names(dn))) # names(NULL) is NULL
           stop("'X' must have named dimnames")
        MARGIN <- match(MARGIN, dnn)
        if (any(is.na(MARGIN)))
            stop("not all elements of 'MARGIN' are names of dimensions")
    }
    
    s.call <- ds [-MARGIN]              # call: with these FUN is called
    s.ans  <- ds [ MARGIN]              # ans:  these stay for the answer
    d.call <- d  [-MARGIN]              # s:    sequence
    d.ans  <- d  [ MARGIN]              # d:    dim
    dn.call<- dn [-MARGIN]              # dn:   dimnames
    dn.ans <- dn [ MARGIN]

    ## do the calls
browser()
    d2 <- prod(d.ans)
    if(d2 == 0L) {
        ## arrays with some 0 extents: return ``empty result'' trying
        ## to use proper mode and dimension:
        ## The following is still a bit `hackish': use non-empty X
        newX <- array(vector(typeof(X), 1L), dim = c(prod(d.call), 1L))
        ans <- FUN(if(length(d.call) < 2L) newX[,1] else
                   array(newX[, 1L], d.call, dn.call), ...)
        return(if(is.null(ans)) ans else if(length(d.ans) < 2L) ans[1L][-1L]
               else array(ans, d.ans, dn.ans))
    }
    ## else
    newX <- aperm(X, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    ans <- vector("list", d2)
    if(length(d.call) < 2L) {# vector
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        for(i in 1L:d2) {
            tmp <- FUN(newX[,i], ...)
            if(!is.null(tmp)) ans[[i]] <- tmp
        }
    } else
       for(i in 1L:d2) {
           tmp <- FUN(array(newX[,i], d.call, dn.call), ...)
           if(!is.null(tmp)) ans[[i]] <- tmp
        }

    ## answer dims and dimnames

    ans.list <- is.recursive(ans[[1L]])
    l.ans <- length(ans[[1L]])

    ans.names <- names(ans[[1L]])
    if(!ans.list)
	ans.list <- any(unlist(lapply(ans, length)) != l.ans)
    if(!ans.list && length(ans.names)) {
        all.same <- vapply(ans, function(x) identical(names(x), ans.names), NA)
        if (!all(all.same)) ans.names <- NULL
    }
    len.a <- if(ans.list) d2 else length(ans <- unlist(ans, recursive = FALSE))
    if(length(MARGIN) == 1L && len.a == d2) {
	names(ans) <- if(length(dn.ans[[1L]])) dn.ans[[1L]] # else NULL
	return(ans)
    }
    if(len.a == d2)
	return(array(ans, d.ans, dn.ans))
    if(len.a && len.a %% d2 == 0L) {
        if(is.null(dn.ans)) dn.ans <- vector(mode="list", length(d.ans))
        dn.ans <- c(list(ans.names), dn.ans)
	return(array(ans, c(len.a %/% d2, d.ans),
		     if(!all(vapply(dn.ans, is.null, NA))) dn.ans))
    }
    return(ans)
}
