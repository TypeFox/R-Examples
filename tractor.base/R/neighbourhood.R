buildStepVectors <- function (width)
{
    smallRes <- width
    largeRes <- smallRes ^ 2
    maxStep <- smallRes %/% 2
    count <- smallRes ^ 3

    directions <- numeric(0)
    for (i in 0:(count-1))
    {
        x <- (i %/% largeRes) - maxStep
        y <- ((i %% largeRes) %/% smallRes) - maxStep
        z <- (i %% smallRes) - maxStep

        directions <- c(directions, x, y, z)
    }

    directionsArray <- array(directions, dim=c(3,count))
    return(directionsArray)
}

#' Image neighbourhoods
#' 
#' This function calculates information about a cuboidal region of an image,
#' with a centre and a fixed voxel width.
#' 
#' @param width An integer voxel width. Must be odd.
#' @param dim An integer giving the dimensionality of the neighbourhood.
#'   Currently must be 3.
#' @param centre A numeric vector giving the centre voxel of the neighbourhood.
#'   Must have exactly \code{dim} elements.
#' @return \code{createNeighbourhoodInfo} returns a list with class
#'   \code{"neighbourhoodInfo"} and elements
#'   \describe{
#'     \item{width}{Copied from the \code{width} argument.}
#'     \item{dim}{Copied from the \code{dim} argument.}
#'     \item{centre}{Copied from the \code{centre} argument.}
#'     \item{vectors}{\code{dim} x \code{width^dim} matrix whose columns give
#'       the locations of each point in the neighbourhood.}
#'     \item{innerProducts}{A square, symmetric matrix of inner products
#'       between every location in the neighbourhood and every other.}
#'   }
#' @author Jon Clayden
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @aliases neighbourhoodInfo
#' @rdname neighbourhoodInfo
#' @export
createNeighbourhoodInfo <- function (width, dim = 3, centre = rep(0,dim))
{
    if (dim != 3)
        report(OL$Error, "Only neighbourhoods in 3 dimensions are supported for now")
    
    vectors <- buildStepVectors(width) + centre
    
    unnormalised <- t(vectors) %*% vectors
    lengths <- sqrt(colSums(vectors^2))
    normalised <- unnormalised / (lengths %o% lengths)
    
    info <- list(width=width, dim=dim, centre=centre, vectors=vectors, innerProducts=normalised)
    class(info) <- "neighbourhoodInfo"
    invisible (info)
}
