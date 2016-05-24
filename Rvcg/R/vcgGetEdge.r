#' Get all edges of a triangular mesh
#'
#' Extract all edges from a mesh and retrieve adjacent faces and vertices
#' @param mesh triangular mesh of class 'mesh3d'
#' @param unique logical: if TRUE each edge is only reported once, if FALSE, all occurences are reported.
#' @return returns a dataframe containing:
#' \item{vert1}{integer indicating the position of the first vertex belonging to this edge}
#' \item{vert2}{integer indicating the position of the second vertex belonging to this edge}
#' \item{facept}{integer pointing to the (or a, if unique = TRUE) face adjacent to the edge}
#' \item{border}{integer indicating if the edge is at the border of the mesh. 0 = no border, 1 = border}
#'
#' @examples
#' require(rgl)
#' data(humface)
#' edges <-vcgGetEdge(humface)
#' \dontrun{
#' ## show first edge
#' lines3d(t(humface$vb[1:3,])[c(edges$vert1[1],edges$vert2[2]),],col=2,lwd=3)
#' shade3d(humface, col=3)
#' ## now find the edge - hint: it is at the neck.
#' }
#' @export
vcgGetEdge <- function(mesh,unique=TRUE)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        mesh <- meshintegrity(mesh,facecheck=TRUE)
        vb <- mesh$vb[1:3,,drop=FALSE]
        it <- mesh$it - 1
        unique <- as.logical(unique)
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        tmp <- .Call("RgetEdge",vb,it,unique)
        edvert <- tmp$edges
        edge <- data.frame(vert1=edvert[,1])
        edge$vert1 <- edvert[,1]
        edge$vert2 <- edvert[,2]
        edge$facept <- tmp$facept
        edge$border <-  tmp$border
        if (!unique)
            edge <- edge[order(edge[,1],edge[,2]),]
        invisible(edge)
    }



#' Get all non-border edges
#'
#' Get all non-border edges and both faces adjacent to them.
#' @param mesh triangular mesh of class 'mesh3d
#' @param silent logical: suppress output of information about number of border edges
#' @return returns a dataframe containing:
#' \item{vert1}{integer indicating the position of the first vertex belonging to this edge}
#' \item{vert2}{integer indicating the position of the second vertex belonging to this edge}
#' \item{border}{integer indicating if the edge is at the border of the mesh. 0 = no border, 1 = border}
#' \item{face1 }{integer pointing to the first face adjacent to the edge}
#' \item{face2 }{integer pointing to the first face adjacent to the edge}
#' @seealso \code{\link{vcgGetEdge}}
#' @examples
#' 
#' data(humface)
#' edges <-vcgNonBorderEdge(humface)
#' ## show first edge (not at the border)
#' \dontrun{
#' require(Morpho)
#' require(rgl)
#' lines3d(t(humface$vb[1:3,])[c(edges$vert1[1],edges$vert2[2]),],col=2,lwd=3)
#' 
#' ## plot barycenters of adjacent faces
#' bary <- barycenter(humface)
#' points3d(bary[c(edges$face1[1],edges$face2[1]),])
#' shade3d(humface, col=3)
#' ## now find the edge - hint: it is at the neck.
#' }
#' @export
vcgNonBorderEdge <- function(mesh, silent=FALSE)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        edges <- vcgGetEdge(mesh,unique=FALSE)
        border <- which(edges$border == 1)
        edgesClean <- edges
        if (length(border) != 0)
            edgesClean <- edgesClean[-border,]
        if (!silent)
            cat(paste("mesh contains ",length(border), "border edges\n"))
        n <- dim(edgesClean)[1]/2
        n2 <- (1:n)*2
        n1 <- n2-1
        out <- edgesClean[n1,]
        out$face1 <- edgesClean$facept[n1]
        out$face2 <- edgesClean$facept[n2]
        out <- out[,-3]

        return(out)
    }
