#' @exportMethod as.SpatialGrid
setGeneric("as.SpatialGrid",
    function(object, nb.cells = 100, cell.size = NULL){ 
        standardGeneric("as.SpatialGrid") 
    }
)

#' Create a spatial grid from an object of class prevR.
#' 
#' This function generates a spatial rectangular grid from the slot \code{boundary} 
#' of an object of class \code{\link[=prevR-class]{prevR}}; function used in particular 
#' by the methods \code{\link[=kde,prevR-method]{kde}}, \code{\link[=krige,prevR-method]{krige}} 
#' and \code{\link[=idw,prevR-method]{idw}}.
#' 
#' @param object object of class \code{\link[=prevR-class]{prevR}}.
#' @param nb.cells number of cells on the longuest side of the studied area 
#'   (unused if \code{cell.size} is defined).
#' @param cell.size size of each cell (in the unit of the projection).
#' 
#' @details This function generates a spatial rectangular grid, each cell being a square 
#' of side \code{cell.size}. If \code{cell.size} is not defined, side of cells will be 
#' calculated as the longuest side of the slot \code{boundary} of \code{object} divided 
#' by \code{nb.cells}.
#' 
#' @return Object of class \code{\link[sp:SpatialGrid-class]{SpatialGrid}}\{\pkg{sp}\}.
#' 
#' @seealso \code{\link[sp]{GridTopology}}\{\pkg{sp}\}, \code{\link[sp]{SpatialGrid-class}}\{\pkg{sp}\}.
#' 
#' @examples 
#' str(as.SpatialGrid(fdhs))
#' str(as.SpatialGrid(fdhs, nb.cells=200))
#'
#' @rdname as.SpatialGrid
#' @aliases as.SpatialGrid-methods as.SpatialGrid,prevR-method as.SpatialGrid
#' @keywords manip spatial

setMethod("as.SpatialGrid","prevR",
  # nb.cells : Un entier qui contient le nombre de cellules sur la plus grande des dimensions (x ou y de clusters)
  #      On deduit facilement la taille d'une cellule donc le nombre de cellules sur la plus petite des dimensions
  # cell.size : la taille d'une cellule. Si cette valeur est fournie nb.cells est ignore
  function (object, nb.cells = 100, cell.size = NULL){
  # calcul de la taille des cellules et donc du nombre de cellules sur chaque coordonnees
  boundary     = slot(object,"boundary")
  if (!attr(boundary,"valid")) {
    clusters  = slot(object,"clusters")
    max.x = max(clusters[["x"]])
    min.x = min(clusters[["x"]])
    max.y = max(clusters[["y"]])
    min.y = min(clusters[["y"]])
  } else {
    bbox  = slot(boundary,"bbox")
    max.x = bbox["x","max"]
    min.x = bbox["x","min"]
    max.y = bbox["y","max"]
    min.y = bbox["y","min"]
  }
  dx = max.x - min.x
  dy = max.y - min.y
  if(is.null(cell.size)){
    dxy = max(c(dx,dy))
    cell.size = dxy/nb.cells
    nb.cells.x = dx%/%cell.size
    nb.cells.y = dy%/%cell.size
  } else {
    nb.cells.x = dx%/%cell.size + 1
    nb.cells.y = dy%/%cell.size + 1
  }

# Creation de la grille sur laquelle aura lieu le lissage
  grid_topo = GridTopology(c(min.x, min.y), c(cell.size, cell.size), c(nb.cells.x, nb.cells.y))
  SG = SpatialGrid(grid_topo, proj4string=object@proj)
  return(SG)
}
)