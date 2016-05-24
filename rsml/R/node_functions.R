
#' Node constructor
#' @param x = x coordinate of the node. Mandatory
#' @param y = y coordinate of the node. Mandatory
#' @param z = z coordinate of the root. Optional
#' @param diameter = diameter of the node. Optional
#' @param orientation = orientation of the node. Optional
#' @param bLength = lenght from the node position in the root from the base of the root. Optional
#' @return the node
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @keywords rsml
#' @export
#' @examples
#' n <- node(1, 1)
node =
  function(x, y, z=0, diameter=0, orientation=0, bLength=0)
  {
    if (!is.numeric(x) || !is.numeric(y) ||
          !all(is.finite(x)) || !all(is.finite(y))
    )
      
      stop("invalid coordinates")
    if (length(x) > 1 || length(y) > 1)
      stop("too big dimension of coordinates")
    nds = list(x = x, y = y, z = z, diameter=diameter, orientation=orientation, bLength=bLength)
    class(nds) = "node"
    nds
  }


#' Print the node
#' @param x object of class node
#' @param ... print options
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @export
#' @examples
#' n <- node(1, 1)
#' print(n)
print.node = 
  function(x, ...)
  {
    obj <- x
    print(paste("x =",format(obj$x)), quote=F)
    print(paste("y = ",format(obj$y)), quote=F)
    print(paste("diameter = ",format(obj$diameter)), quote=F)
    print(paste("orientation = ",format(obj$orientation)), quote=F)
    print(paste("distance from base = ",format(obj$bLength)), quote=F)
  }
