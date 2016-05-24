#' For mapping global data.
#' 
#' Enables mapping of country level and gridded user datasets by facilitating
#' joining to modern world maps and offering visualisation options. Country
#' borders are derived from Natural Earth data v 1.4.0.
#' 
#' Country Level Data can be joined to a map using
#' \code{\link{joinCountryData2Map}}, then mapped using
#' \code{\link{mapCountryData}}. These functions can cope with a range of
#' country names and country codes.
#' 
#' Country boundaries are derived from version 1.4.0 of Natural Earth data as
#' described in \code{\link{countriesCoarse}}. Higher resolution boundaries are
#' provided in a companion package rworldxtra.
#' 
#' More generic functions allow the user to provide their own polygon map using
#' \code{\link{joinData2Map}} and \code{\link{mapPolys}}.
#' 
#' Bubble, bar and pie charts can be added to maps using
#' \code{\link{mapBubbles}}, \code{\link{mapBars}} and \code{\link{mapPies}}.
#' 
#' Try the new method \code{\link{barplotCountryData}} for producing a ranked
#' bar plot of country data with country names that can provide a useful
#' companion to maps.
#' 
#' Options are provided for categorising data, colouring maps and symbols, and
#' adding legends.
#' 
#' Gridded data can be mapped using \code{\link{mapGriddedData}}, but the
#' raster package is much more comprehensive.
#' 
#' Type vignette('rworldmap') to access a short document showing a few examples
#' of the main rworldmap functions to get you started.
#' 
#' @name rworldmap-package
#' @aliases rworldmap-package rworldmap
#' @docType package
#' @author Andy South
#' 
#' with contributions from Joe Scutt-Phillips, Barry Rowlingson, Roger Bivand
#' and Pru Foster
#' 
#' Maintainer: <southandy@@gmail.com>
#' @references Stable version :
#' http://cran.r-project.org/web/packages/rworldmap 
#' \cr Development version :https://github.com/AndySouth/rworldmap
#' 
#' Discussion group : http://groups.google.com/group/rworldmap
#' @keywords package
#' @import sp grDevices graphics
#' @importFrom methods is
#' @importFrom stats aggregate quantile
#' @importFrom utils data str

NULL



