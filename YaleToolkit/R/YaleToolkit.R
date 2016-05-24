#' Data exploration tools from the Department of Statistics at Yale University
#' 
#' This collection of data exploration tools was developed
#' at Yale University for the graphical exploration of complex
#' multivariate data. The main functions provided are \code{barcode()},
#' \code{gpairs()}, \code{whatis()}, and \code{sparkmat()}, although \code{barcode()}
#' and \code{gpairs()} are now provided by packages of the same names, respectively.
#'
#' The package also includes several data sets.  For more information,
#' please see the help files for \code{nasa} and \code{YaleEnergy}.
#' Please get in touch with us if you note any problems.
#' 
#' @references
#' \itemize{
#'   \item Chambers, J.M., Cleveland, W.S., Kleiner, B., and Tukey, P.A. (1983),
#'     \emph{Graphical Methods for Data Analysis}, Belmont, CA: Wadsworth.
#'   \item Friendly, M. (2002) 'Corrgrams: Exploratory displays for correlation
#'     matrices' \code{American Statistician} 56(4), 316--324.
#'   \item Tufte, Edward R. (2006) \emph{Beautiful Evidence} The Graphics
#'     Press, Cheshire, Connecticut.
#'     See \url{http://www.edwardtufte.com} for this and other references.
#' }
#' 
#' @author John W. Emerson, Walton Green
#' @docType package
#' @name YaleToolkit
#' @export sparkline
#' @export sparklines
#' @export sparkmat
#' @importFrom foreach foreach
#' @importFrom foreach '%do%'
#' @importFrom iterators iread.table
#' @importFrom iterators ireadLines
#' @import grid
NULL