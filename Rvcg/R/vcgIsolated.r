#' Remove isolated pieces from a surface mesh or split into connected components
#' 
#' Remove isolated pieces from a surface mesh, selected by a
#' minimum amount of faces or of a diameter below a given threshold.
#' Also the option only to keep the largest piece can be selected or to split a mesh into connected components.
#' 
#' 
#' @param mesh triangular mesh of class "mesh3d".
#' @param facenum integer: all connected pieces with less components are
#' removed. If not specified or 0 and diameter is NULL, then only the component
#' with the most faces is kept. 
#' @param diameter numeric: all connected pieces smaller diameter are removed
#' removed. \code{diameter = 0} removes all component but the largest ones. This option overrides the option \code{facenum}.
#' @param split logical: if TRUE, a list with all connected components of the mesh will be returned.
#' @param silent logical, if TRUE no console output is issued.
#' @return returns the reduced mesh.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(humface)
#' cleanface <- vcgIsolated(humface)
#' 
#' 
#' @export vcgIsolated
vcgIsolated <- function(mesh,facenum=NULL,diameter=NULL,split=FALSE,silent=FALSE) {
    if (!inherits(mesh,"mesh3d"))
        stop("argument 'mesh' needs to be object of class 'mesh3d'")
    if (is.null(facenum))
        facenum <- 0
    if (!is.null(diameter))
        facenum <- -1

    if (is.null(diameter))
        diameter <- 0
    storage.mode(facenum) <- "integer"
    storage.mode(diameter) <- "double"

    mesh <- meshintegrity(mesh,facecheck=TRUE)
    vb <- mesh$vb[1:3,,drop=FALSE]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    tmp <- .Call("Risolated", vb, it, diameter, facenum,silent,split)
    if (!split) {
        ## handle vertex color
        if (!is.null(mesh$material$color)) {
            if (length(tmp$remvert)) {
                colframe <- data.frame(it=1:ncol(mesh$vb))
                colframe$rgb <- rep("#FFFFFF",ncol(mesh$vb))
                colframe$it <- 1:ncol(mesh$vb)
                remvert <- tmp$remvert
                tmp1 <- data.frame(it=as.vector(mesh$it))
                tmp1$rgb <- as.vector(mesh$material$color)
                tmp1 <- unique(tmp1)
                tmp1 <- tmp1[order(tmp1$it),]
                colframe$rgb[tmp1$it] <- tmp1$rgb
                colvec <- colframe$rgb[!as.logical(remvert)]
                colfun <- function(x) {
                    x <- colvec[x]
                    return(x)
                }
                tmp$material$color <- matrix(colfun(tmp$it),dim(tmp$it))
            } else {
                tmp$material$color <- mesh$material$color
            }
        }
    }
    
    invisible(tmp)
}
