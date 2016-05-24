#' @rdname vector
#' @export
resolveVector <- function (len, ...)
{
    vector <- c(...)
    if (is.numeric(vector) && (length(vector) == len))
        return (vector)
    else
        return (NULL)
}

#' Miscellaneous vector functions
#' 
#' These functions provide the (Euclidean) length of a vector, the vector cross
#' product or angle between two vectors.
#' 
#' @param vector,v1,v2 Numeric vectors of any length.
#' @param a,b Numeric 3-vectors.
#' @param len The expected length of the vector.
#' @param \dots Elements of the vector, to be concatenated together.
#' @return For \code{vectorLength}, the Euclidean norm or length of the
#'   specified vector, given by \code{sqrt(sum(vector^2))}. For
#'   \code{vectorCrossProduct}, the vector cross product of the two specified
#'   vectors; and for \code{angleBetweenVectors}, the angle (in radians)
#'   between the two specified vectors. The \code{resolveVector} function
#'   concatenates the values given in \code{\dots{}}, and if the result is a
#'   vector of length \code{len} then it is returned. If not, \code{NULL} is
#'   returned.
#' @author Jon Clayden
#' @seealso \code{\link{crossprod}} for the matrix cross product.
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @rdname vector
#' @export
vectorLength <- function (vector)
{
    return (sqrt(sum(vector^2)))
}

#' @rdname vector
#' @export
vectorCrossProduct <- function (a, b)
{
    if (length(a) != 3 || length(b) != 3)
        report(OL$Error, "Cross product is currently only defined for 3-vectors")
    
    # Ref: http://mathworld.wolfram.com/CrossProduct.html
    return (c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]))
}

#' @rdname vector
#' @export
angleBetweenVectors <- function (v1, v2)
{
    if (is.na(v1) || is.na(v2))
        return (NA)
    else
    {
        if (identical(v1,v2))
            cosine <- 1
        else
            cosine <- (v1 %*% v2) / (vectorLength(v1) * vectorLength(v2))
        return (acos(cosine))
    }
}

matrixToVectorLocs <- function (matrixLocs, dims)
{
    nDims <- length(dims)
    storage.mode(matrixLocs) <- "integer"
    matrixLocs <- promote(matrixLocs, byrow=TRUE)
    jumps <- c(1, cumprod(dims))
    return (rowSums((matrixLocs - 1) * rep(jumps[1:nDims],each=nrow(matrixLocs))) + 1)
}

vectorToMatrixLocs <- function (vectorLocs, dims)
{
    return (arrayInd(vectorLocs, dims))
}
