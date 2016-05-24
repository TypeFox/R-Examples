mtapply <-
function (X, INDEX, FUN = NULL, simplify = TRUE)
{
    FUN <- if (!is.null(FUN))
        match.fun(FUN)
    if (!is.list(INDEX))
        INDEX <- list(INDEX)
    nI <- length(INDEX)
    namelist <- vector("list", nI)
    names(namelist) <- names(INDEX)
    extent <- integer(nI)
	
	nx <- ifelse(is.list(X), 
					ifelse(is.data.frame(X[[1]]), dim(X[[1]])[1], length(X[[1]])), 
					length(X))
	
    one <- as.integer(1)
    group <- rep.int(one, nx)
    ngroup <- one
    for (i in seq(INDEX)) {
        index <- as.factor(INDEX[[i]])
        if (length(index) != nx)
            stop("arguments must have same length")
        namelist[[i]] <- levels(index)
        extent[i] <- nlevels(index)
        group <- group + ngroup * (as.integer(index) - one)
        ngroup <- ngroup * nlevels(index)
    }
    if (is.null(FUN))
        return(group)
     if (!is.list(X)) {
    ans <- lapply(split(X, group), FUN)
    index <- as.numeric(names(ans))
     }
     else {
    myargs<-vector("list",length(X)+1)
     for (i in 1:length(X)) myargs[[i+1]]<-split(X[[i]],group)
     myargs[[1]]<-FUN
     ans<-do.call(mapply,myargs)
     ansx <- lapply(myargs[[2]],length)
     index <- as.numeric(names(ansx))
     }
    if (simplify && all(unlist(lapply(ans,length)) == 1)) {
        ansmat <- array(dim = extent, dimnames = namelist)
        if (is.list(ans)) ans <- unlist(ans, recursive = FALSE)
    }
    else {
        ansmat <- array(vector("list", prod(extent)), dim = extent,
            dimnames = namelist)
    }
    names(ans) <- NULL
    ansmat[index] <- ans
    ansmat
}

