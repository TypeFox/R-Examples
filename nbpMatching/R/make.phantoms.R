#'Add Phantom Rows and Columns
#'
#'The make.phantoms function will take an \eqn{N}x\eqn{N} matrix and add
#'\eqn{NP} phantom elements, thus creating a matrix with \eqn{N+NP}x\eqn{N+NP}
#'dimensions.
#'
#'This function is internal to the \code{\link{gendistance}} function, but may be
#'useful in manufacturing personalized distance matrices.  Phantoms are fake
#'elements that perfectly match all elements.  They can be used to discard a
#'certain number of elements.
#'
#'@aliases make.phantoms make.phantoms,ANY,missing-method
#'make.phantoms,data.frame,numeric-method make.phantoms,matrix,numeric-method
#'@param x A matrix or data.frame object, with \eqn{N}x\eqn{N} dimensions.
#'@param nphantoms An integer, providing the number of phantom elements to add.
#'@param name A character string, indicating the name attribute for new
#'elements.  Defaults to "phantom".
#'@param maxval An integer value, the default value to give the pairs of
#'phantoms (indeces [N+1:N+NP, N+1:N+NP]), assumed to be a maximum distance.
#'Defaults to Inf.
#'@param \dots Additional arguments, not used at this time.
#'@return a matrix or data.frame object
#'@exportMethod make.phantoms
#'@author Cole Beck
#'@seealso \code{\link{gendistance}} \code{\link{distancematrix}}
#'@examples
#'
#'# 5x5 distance matrix
#'dist.mat <- matrix(c(0,5,10,15,20,5,0,15,25,35,10,15,0,25,40,15,25,25,0,15,20,35,40,15,0), nrow=5)
#'# add one phantom element
#'dm.ph <- make.phantoms(dist.mat, 1)
#'# create distancematrix object
#'distancematrix(dm.ph)
#'# add three phantoms
#'make.phantoms(dist.mat, 3)
#'

setGeneric("make.phantoms", function(x, nphantoms, name="phantom", maxval=Inf, ...) standardGeneric("make.phantoms"))
setMethod("make.phantoms", signature(x="matrix", nphantoms="numeric"), function(x, nphantoms, name="phantom", maxval=Inf, ...) {
    if(nphantoms < 1) return(x)
    if(missing(name)) {
        name <- "phantom"
    } else if(!is.character(name)) {
        stop("name argument is not character")
    }
    if(missing(maxval)) {
        maxval <- Inf
    } else if(!is.numeric(maxval)) {
        stop("maxval argument is not numeric")
    }
    nr <- nrow(x)
    # preserve rownames
    mynames <- rownames(x)
    if(is.null(mynames)) {
        mynames <- seq(nr)
    }
    # phantom index
    p.index <- seq(from=nr+1, length.out=nphantoms)
    newvals <- rep(0, nphantoms)
    # add phantom columns
    m <- do.call("cbind", c(list(x), newvals))
    # add phantom rows
    m <- do.call("rbind", c(list(m), newvals))
    # distance between phantoms should be maxval
    m[p.index, p.index] <- maxval
    # create row names
    mynames <- c(mynames, sprintf("%s%s", name, p.index))
    dimnames(m) <- list(mynames, mynames)
    m
})
# x is data.frame instead of matrix
setMethod("make.phantoms", signature(x="data.frame", nphantoms="numeric"), function(x, nphantoms, name, maxval, ...) {
    # convert to matrix and back
    as.data.frame(make.phantoms(as.matrix(x), nphantoms, name, maxval, ...))
})
# don't do anything when nphantoms is missing
setMethod("make.phantoms", signature(nphantoms="missing"), function(x, nphantoms, name, maxval, ...) {
    return(x)
})
