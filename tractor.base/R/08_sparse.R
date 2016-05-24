#' The SparseArray class
#' 
#' This class represents an array with any number of dimensions, in which a
#' significant proportion of entries are zero. The coordinates of nonzero
#' entries are stored along with their values, with all remaining entries
#' assumed to be zero. Methods are provided to index into the array in the
#' standard way, using matrix or vector indices; and for coercing between
#' \code{SparseArray} objects and standard (dense) arrays.
#' 
#' @field data Vector of nonzero data values
#' @field coords Integer matrix of nonzero \code{data} locations, one per row
#' @field dims Integer vector of dimensions
#' 
#' @export
SparseArray <- setRefClass("SparseArray", contains="SerialisableObject", fields=list(data="ANY",coords="matrix",dims="integer"), methods=list(
    initialize = function (...)
    {
        object <- initFields(...)
        storage.mode(object$coords) <- "integer"
        return (object)
    },
    
    aperm = function (perm)
    {
        "Permute the dimensions of the array"
        .self$coords <- .self$coords[,perm]
        .self$dims <- .self$dims[perm]
    },
    
    apply = function (margin, fun, ...)
    {
        "Apply a function to margins of the array"
        fun <- match.fun(fun)
        dimsToKeep <- .self$dims[margin]
        dimsToLose <- .self$dims[-margin]
        iterations <- prod(dimsToKeep)
        
        reshapedArray <- .self$copy()
        reshapedArray$aperm(c(setdiff(seq_len(.self$getDimensionality()),margin), margin))
        reshapedArray$setDimensions(c(prod(dimsToLose), iterations))
        
        # Assume for now that each application of the function will result in a vector of the same length
        tempData <- as.array(reshapedArray[,1])
        dim(tempData) <- dimsToLose
        tempResult <- fun(tempData, ...)
        finalArray <- array(NA, dim=c(length(tempResult),iterations))
        for (i in 1:iterations)
        {
            tempData <- as.array(reshapedArray[,i])
            dim(tempData) <- dimsToLose
            finalArray[,i] <- fun(tempData, ...)
        }
        
        dim(finalArray) <- c(dim(finalArray)[1], dimsToKeep)
        finalArray <- drop(finalArray)
        if (length(dim(finalArray)) == 1)
            dim(finalArray) <- NULL
        
        return (finalArray)
    },
    
    flip = function (dimsToFlip)
    {
        "Flip the array along one or more directions"
        for (i in dimsToFlip)
            .self$coords[,i] <- .self$dims[i] - .self$coords[,i] + 1
    },
    
    getCoordinates = function () { return (.self$coords) },
    
    getData = function () { return (.self$data) },
    
    getDimensionality = function () { return (length(.self$dims)) },
    
    getDimensions = function () { return (.self$dims) },
    
    setCoordinatesAndData = function (newCoords, newData)
    {
        "Update the nonzero locations and data values in the array"
        .self$coords <- newCoords
        .self$data <- newData
    },
    
    setDimensions = function (newDims)
    {
        "Change the dimensions of the image"
        if (prod(.self$dims) != prod(newDims))
            report(OL$Error, "New dimensions are incompatible with this SparseArray object")
        .self$coords <- vectorToMatrixLocs(matrixToVectorLocs(.self$coords,.self$dims), newDims)
        .self$dims <- as.integer(newDims)
    }
))

.evaluateIndices <- function (i, j, ...)
{
    # Find missing arguments and replace with NULL (code borrowed from the "slam" package)
    args <- substitute(list(i, j, ...), parent.frame())
    argsEmpty <- sapply(args, function(a) identical(as.character(a), ""))
    args[argsEmpty] <- list(NULL)
    return (eval(args, envir=parent.frame(2)))
}

#' Indexing methods
#' 
#' Indexing methods for \code{\link{SparseArray}} and \code{\link{MriImage}}
#' objects. For the latter class, arguments are passed to the equivalents for
#' \code{array} or \code{\link{SparseArray}}. For \code{\link{SparseArray}},
#' indexing may be blank, or by numeric vector or matrix.
#' 
#' @param x An object of the appropriate type.
#' @param i,j,\dots Indexing objects.
#' @param drop Scalar value: should unitary dimensions be dropped?
#' @param value New value(s) for replacement forms.
#' @return A vector, array or \code{\link{SparseArray}}.
#' 
#' @author Jon Clayden
#' @rdname index
#' @export
setMethod("[", "SparseArray", function (x, i, j, ..., drop = TRUE) {
    # This implementation owes a lot to the equivalent in the "slam" package (credit: Kurt Hornik, David Meyer and Christian Buchta)
    nArgs <- nargs() - as.integer(!missing(drop))
    
    nDims <- x$getDimensionality()
    dims <- x$getDimensions()
    data <- x$getData()
    coords <- x$getCoordinates()
    
    if (nArgs < 2)
        return (data)
    else if (nArgs == 2)
    {
        if (is.list(i))
            args <- i
        else
        {
            index <- i
            if (is.logical(index) || is.character(index))
                report(OL$Error, "Logical and character indexing are not yet supported")
            else if (is.matrix(index))
            {
                # Matrix indexing, one row per point: convert to vector and drop through
                if (ncol(index) != nDims)
                    report(OL$Error, "Number of dimensions given does not match image")
                index <- matrixToVectorLocs(index, dims)
            }
        
            # Vector indexing, one number per point
            if (any(index <= 0))
                report(OL$Error, "Zero and negative indices are not yet supported")
        
            returnValue <- vector(mode=storage.mode(data), length=length(index))
            dataLocs <- match(index, matrixToVectorLocs(coords,dims), 0L)
            returnValue[dataLocs > 0] <- data[dataLocs]
        
            return (returnValue)
        }
    }
    else if (nArgs != (nDims + 1))
        report(OL$Error, "Number of dimensions given does not match image")
    else
        args <- .evaluateIndices(i, j, ...)
    
    dataToKeep <- rep.int(TRUE, length(data))
    finalDims <- dims
    finalCoords <- coords
    
    for (currentDim in 1:nDims)
    {
        currentDimIndex <- args[[currentDim]]
        if (is.null(currentDimIndex))
            next
        else if (!is.numeric(currentDimIndex))
            report(OL$Error, "Only numeric indices are currently supported")
        else if (any(currentDimIndex <= 0))
            report(OL$Error, "Zero and negative indices are not yet supported")
        else
        {
            finalDims[currentDim] <- length(currentDimIndex)
            dataLocs <- match(coords[,currentDim], currentDimIndex, 0L)
            finalCoords[dataLocs > 0, currentDim] <- seq_along(currentDimIndex)[dataLocs]
            dataToKeep <- dataToKeep & (dataLocs > 0)
        }
    }
    
    if (drop)
        dimsToKeep <- which(finalDims > 1)
    else
        dimsToKeep <- 1:nDims
    
    if (length(dimsToKeep) < 2)
    {
        # Only one dimension is nonunitary, so the index into the final
        # vector is the maximum of the indices into each dimension for
        # each nonzero data value - the other indices will all be 1
        returnValue <- vector(mode=storage.mode(data), length=prod(finalDims))
        returnValue[apply(finalCoords[dataToKeep,,drop=FALSE],1,max)] <- data[dataToKeep]
    }
    else
        returnValue <- SparseArray$new(data=data[dataToKeep], coords=finalCoords[dataToKeep,dimsToKeep,drop=FALSE], dims=finalDims[dimsToKeep])
    
    return (returnValue)
})

#' @rdname index
#' @export
setReplaceMethod("[", "SparseArray", function (x, i, j, ..., value) {
    nArgs <- nargs() - 1
    
    nDims <- x$getDimensionality()
    dims <- x$getDimensions()
    data <- x$getData()
    coords <- x$getCoordinates()
    
    if (nArgs < 2)
        index <- matrixToVectorLocs(coords, dims)
    else if (nArgs == 2 && !is.list(i))
    {
        index <- i
        if (is.logical(index) || is.character(index))
            report(OL$Error, "Logical and character indexing are not yet supported")
        else if (is.matrix(index))
        {
            # Matrix indexing, one row per point: convert to vector and drop through
            if (ncol(index) != nDims)
                report(OL$Error, "Number of dimensions given does not match image")
            index <- matrixToVectorLocs(index, dims)
        }
        
        # Vector indexing, one number per point
        if (any(index <= 0))
            report(OL$Error, "Zero and negative indices are not yet supported")
    }
    else if (!is.list(i) && nArgs != (nDims + 1))
        report(OL$Error, "Number of dimensions given does not match image")
    else
    {
        if (is.list(i))
            args <- i
        else
            args <- .evaluateIndices(i, j, ...)
        
        args <- lapply(seq_len(nDims), function(i) {
            if (is.null(args[i]))
                1:dims[i]
            else
                args[i]
        })
        
        index <- as.matrix(expand.grid(args))
        storage.mode(index) <- "integer"
        index <- matrixToVectorLocs(index, dims)
    }
    
    if (length(index) %% length(value) != 0)
        report(OL$Error, "Number of items to replace is not a multiple of replacement length")
    if (length(index) != length(value))
        value <- rep(value, length(index) %/% length(value))
    
    dataLocs <- match(index, matrixToVectorLocs(coords,dims), 0L)
    present <- (dataLocs > 0)
    zero <- (value == 0)
    
    data[which(present & !zero)] <- value[which(present & !zero)]
    if (any(present & zero))
    {
        coords <- coords[-which(present & zero),]
        data <- data[-which(present & zero)]
    }
    if (any(!present & !zero))
    {
        coords <- rbind(coords, arrayInd(index[which(!present & !zero)],dims))
        data <- c(data, value[which(!present & !zero)])
    }
    
    x$setCoordinatesAndData(coords, data)
    return (x)
})

setMethod("Ops", signature(e1="SparseArray",e2="ANY"), function (e1, e2) {
    report(OL$Error, "No operator method is defined for a \"SparseArray\" and \"#{class(e2)[1]}\"")
})

setMethod("Ops", signature(e1="SparseArray",e2="array"), function (e1, e2) {
    if (!all(callGeneric(0, unique(e2)) == 0))
        return (callGeneric(as(e1,"array"), e2))
    else
    {
        if (!equivalent(e1$getDimensions(), dim(e2)))
            report(OL$Error, "Sparse and dense array dimensions don't match")
    
        newData <- callGeneric(e1$getData(), e2[e1$getCoordinates()])
        zero <- (newData == 0)
        return (newSparseArrayWithData(newData[!zero], e1$getCoordinates()[!zero,,drop=FALSE], e1$getDimensions()))
    }
})

setMethod("Ops", signature(e1="SparseArray",e2="numeric"), function (e1, e2) {
    if (callGeneric(0, e2) != 0)
        return (callGeneric(as(e1,"array"), e2))
    else
    {        
        newData <- callGeneric(e1$getData(), e2)
        zero <- (newData == 0)
        return (newSparseArrayWithData(newData[!zero], e1$getCoordinates()[!zero,,drop=FALSE], e1$getDimensions()))
    }
})

setMethod("Summary", signature(x="SparseArray"), function(x, ..., na.rm = FALSE) {
    callGeneric(0, x$getData(), ..., na.rm=na.rm)
})

setAs("array", "SparseArray", function (from) {
    coordinates <- which(!is.na(from) & from != 0, arr.ind=TRUE)
    object <- SparseArray$new(data=from[coordinates], coords=coordinates, dims=dim(from))
    return (object)
})

setAs("SparseArray", "array", function (from) {
    data <- array(vector(mode=storage.mode(from$getData()),length=1), dim=from$getDimensions())
    data[from$getCoordinates()] <- from$getData()
    return (data)
})

setAs("SparseArray", "logical", function (from) {
    data <- array(FALSE, dim=from$getDimensions())
    data[from$getCoordinates()] <- TRUE
    return (as.vector(data))
})

#' @export
as.array.SparseArray <- function (x, ...)
{
    as(x, "array")
}

#' @export
as.vector.SparseArray <- function (x, mode = "any")
{
    as.vector(as(x,"array"), mode=mode)
}

#' @export
dim.SparseArray <- function (x)
{
    x$getDimensions()
}

#' @export
"dim<-.SparseArray" <- function (x, value)
{
    x$setDimensions(value)
    return (x)
}

#' @export
aperm.SparseArray <- function (a, perm, ...)
{
    newObject <- a$copy()$aperm(perm)
    return (newObject)
}

#' Create a SparseArray object
#' 
#' This function creates a \code{\linkS4class{SparseArray}} object from its
#' constituent parts.
#' 
#' @param data A vector of (nonzero) array elements.
#' @param coordinates A matrix with as many rows as \code{data} has elements,
#'   containing the coordinates of each nonzero element in the array.
#' @param dims The dimensions of the array.
#' @return A \code{\linkS4class{SparseArray}} object.
#' @author Jon Clayden
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export newSparseArrayWithData
newSparseArrayWithData <- function (data, coordinates, dims)
{
    invisible (SparseArray$new(data=data, coords=coordinates, dims=dims))
}
