#' Clean triangular surface meshes
#'
#' Apply several cleaning algorithms to surface meshes
#' @param mesh triangular mesh of class 'mesh3d'
#' @param sel integer vector selecting cleaning type (see "details"), 
#' @param tol numeric value determining Vertex Displacement Ratio used for splitting non-manifold vertices.
#' @param silent logical, if TRUE no console output is issued.
#' @param iterate logical: if TRUE, vcgClean is repeatedly run until nothing more is to be cleaned (see details).
#' @details the vector sel determines which operations are performed in which order. E.g. removing degenerate faces may generate unreferenced vertices, thus the ordering of cleaning operations is important, multiple calls are possible (sel=c(1,3,1) will remove unreferenced vertices twice).
#' available options are:
#'
#' \itemize{
#' \item{0 = only duplicated vertices and faces are removed}
#' \item{1 = remove unreferenced vertices}
#' \item{2 = Remove non-manifold Faces}
#' \item{3 = Remove degenerate faces}
#' \item{4 = Remove non-manifold vertices}
#' \item{5 = Split non-manifold vertices by threshold}
#' \item{6 = merge close vertices (radius=\code{tol})}
#' \item{7 = coherently orient faces}
#' }
#' @return cleaned mesh with an additional entry
#' \item{remvert}{vector of length = number of vertices before cleaning. Entries = 1 indicate that this vertex was removed; 0 otherwise.}
#' @examples
#' data(humface)
#' cleanface <- humface
#' ##add duplicated faces
#' cleanface$it <- cbind(cleanface$it, cleanface$it[,1:100])
#' ## add duplicated vertices
#' cleanface$vb <- cbind(cleanface$vb,cleanface$vb[,1:100])
#' ## ad unreferenced vertices
#' cleanface$vb <- cbind(cleanface$vb,rbind(matrix(rnorm(18),3,6),1))
#' cleanface <- vcgClean(cleanface, sel=1)
#' @export vcgClean
vcgClean <- function(mesh, sel = 0,tol=0,silent=FALSE,iterate=FALSE) {
    removed <- 10
    
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
    while (removed > 0) {
        mesh <- meshintegrity(mesh)
        vb <- mesh$vb[1:3,,drop=FALSE]
        it <- mesh$it - 1
        tol <- tol[1]
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        sel <- as.vector(sel)
        tmp <- .Call("Rclean", vb, it, sel, tol,silent)
        removed <- sum(tmp$remvert)
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
        mesh <- tmp
        if (!iterate)
            removed <- 0
    }
    
    return(meshintegrity(mesh))
    
}
