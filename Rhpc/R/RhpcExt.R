Rhpc_EvalQ<-function(cl, expr, envir=.GlobalEnv)
    Rhpc_worker_call(cl, eval, substitute(expr), envir=envir)


Rhpc_apply<- function(cl = NULL, X, MARGIN, FUN, ...)
{
    ## from parApply in parallel package
    ## rewrite later

    FUN <- match.fun(FUN) # should this be done on worker?

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
        if (anyNA(MARGIN))
            stop("not all elements of 'MARGIN' are names of dimensions")
    }
    s.call <- ds[-MARGIN]
    s.ans  <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans  <- d[MARGIN]
    dn.call <- dn[-MARGIN]
    dn.ans <- dn[MARGIN]
    ## dimnames(X) <- NULL

    ## do the calls

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
    arglist <- if(length(d.call) < 2L) {# vector
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        lapply(seq_len(d2), function(i) newX[,i])
    } else
        lapply(seq_len(d2), function(i) array(newX[,i], d.call, dn.call))
    ans <- Rhpc_lapply(cl = cl, X = arglist, FUN = FUN, ...)

    ## answer dims and dimnames

    ans.list <- is.recursive(ans[[1L]])
    l.ans <- length(ans[[1L]])

    ans.names <- names(ans[[1L]])
    if(!ans.list)
        ans.list <- any(lengths(ans) != l.ans)
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

Rhpc_sapply<-function (cl = NULL, X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) 
{
    ## rewrite later

    FUN <- match.fun(FUN)
    answer <- Rhpc_lapply(cl, X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!identical(simplify, FALSE) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

Rhpc_sapplyLB<-function (cl = NULL, X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) 
{
    ## rewrite later

    FUN <- match.fun(FUN)
    answer <- Rhpc_lapplyLB(cl, X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!identical(simplify, FALSE) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}
