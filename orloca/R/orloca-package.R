#' The package deals with Operations Research LOCational Analysis models
#' 
#' This version of the package deals with the min-sum location problem, also known as Fermat--Weber problem.
#' 
#' The min-sum location problem look for a point such that the weighted sum of the distances to the demand points are minimized.
#'
#' @aliases orloca-package
#' @docType package
#' @name orloca-package
#' @details
#' \preformatted{
#' 
#' Package:   orloca
#' 
#' Type:      Package
#' 
#' Version:   4.2
#' 
#' Date:      2014-06-02
#' 
#' License:   GPL (>= 3)
#' }
#'
#' The package provides a class (\code{loca.p}) that represents a location problem with a finite set of demand points over the plane.
#' Also, it is possible to plot the points and the objective function.
#' Such objective function is the total weighted distances travelled by all the customers to the service.
#'
#'
#' Non-planar location problems could be handle in future versions of the package.
#'
#' 
#' For a demo, load the package with \code{library(orloca)}, and use \code{demo(orloca)}.
#'
#' 
#' The package is ready for internationalization. The authors ask for translated version of the .mo file to include in the package.
#'
#' @author Fernando Fernandez-Palacin <fernando.fernandez@@uca.es> and Manuel Munoz-Marquez <manuel.munoz@@uca.es>
#' 
#' Mantainer: Manuel Munoz-Marquez <manuel.munoz@@uca.es>
#' @references
#' [1] Love, R. F., Morris, J. G., Wesolowsky, G. O. \emph{Facilities Location: Chapter 2: Introduction to Single-Facility Location}, 1988, North-Holland
#'
#' [2] \url{http://knuth.uca.es/orloca}
#' @keywords package optimize
#' @seealso
#' Para la version en espanol, instale el paquete orloca.es y consulte la ayuda sobre \code{\link[orloca.es]{orloca.es-package}}.
#' (For the spanish version, install the package orloca.es and see the help about \code{\link[orloca.es]{orloca.es-package}}).
#' @examples
#' # A new unweighted loca.p object
#' o <- loca.p(x = c(-1, 1, 1, -1), y = c(-1, -1, 1, 1))
#' 
#' # Compute the sum of distances to point (3, 4)
#' zsum(o, 3, 4)
#' 
#' # Compute the sum of distances to point (3, 4) using lp norm
#' zsum(o, 3, 4, lp=2.5)
#'
#' # Solve the optimization problem
#' zsummin(o)
#' # Contour plot
#' contour(o)
#'
#' # Make a demo of the package
#' demo(orloca)
#' 
#' @import methods
#'
#' @export as.loca.p
#' @export as.loca.p.matrix
#' @export as.matrix.loca.p
#' @export contour.loca.p
#' @export loca.p
#' @export persp.loca.p
#' @export plot.loca.p
#' @export zsum
#' @export zsumgra
#' @export zsummin
#'
NULL
