#' Convert a surface in xyz to a data frame.
#' 
#' Several functions (for example \code{\link[GenKern]{KernSur}}\{\pkg{GenKern}\}) 
#' return a surface as a list "xyz" composed of three elements: vector of ordinates in the x dimension, 
#' vector of ordinates in the y dimension and a matrix with the values of the surface in x and y. 
#' This function transforms a list "xyz" into a data frame.
#' 
#' @param xyz a list with 3 elements: a vector with x-coordinates, a vector with y-coordinates and 
#'   and matrix with value for each point of coordinates x[i],y[j].
#' @param xcol x index.
#' @param ycol y index.
#' @param zcol z index.
#' 
#' @note \code{xyz} could be a list like \code{x,y,z1,z2,z3}. If so, \code{zcol} should be equal 
#' to \code{c("z1","z2","z3")} or \code{c(3,4,5)}.
#' 
#' @return A \code{data.frame}.
#' 
#' @examples 
#'   x <- c(2,4,6,8,10)
#'   y <- x 
#'   op <- GenKern::KernSur(x,y, xgridsize=50, ygridsize=50,
#'                 correlation=0, 
#'                 xbandwidth=1, ybandwidth=1,
#'                 range.x=c(0,13), range.y=c(0,13)
#'   )
#'   str(op)
#'   
#'   op.df <- xyz2dataframe(op)
#'   str(op.df)
#'
#' @keywords manip spatial
#' @export

# Transforme une liste de 3 elements (x, y et z) comme le resultat de KernSur
# en un data.frame
# note : z peut etre une liste de matrices
xyz2dataframe <- function(xyz, xcol=1, ycol=2, zcol=3) {
    nx <- length(xyz[[xcol]])
    ny <- length(xyz[[ycol]])
    for (i in 1:ny) {
        z <- xyz[zcol]
        for (k in 1:length(zcol)) z[[k]] <- z[[k]][,i]
        y <- rep(xyz[[ycol]][i],nx)
        x <- xyz[xcol]
        temp <- data.frame(x,y,z)
        names(temp) <- c(names(xyz[xcol]),names(xyz[ycol]),names(xyz[zcol]))
        if (i==1)
            res <- temp
        else
            res <- rbind(res,temp)
    }
    return(res)
}