#' @name aoristic-package
#' @docType package
#' @title Creating a kml file with aoristic graph: Aoristic Analysis with spatial kml data
#' @description Using crime incident data with lat/lon, DateTimeFrom, and DateTimeTo, a total of three (3) aoristic graphs can be created: 1) density and contour; 2) grid count; and 3) shapefile boundary.
#' @details 
#' Google Earth needs to be installed in order to view a kml output with aoristic graphs.  
#' aoristic.grid produces a kml file of aoristic graphs based on a 5*5 grid.
#' aoristic.shp produces a kml file of aoristic graphs based on input shapefile boundary
#' aoristic.density produces a kml file of aoristic graphs with kernel density
#' @author George Kikuchi \email{gkikuchi@@csufresno.edu}
#' @references Ratcliffe, J. H. (2002). Aoristic Signatures and the Spatio-Temporal Analysis of High Volume Crime Patterns. Journal of Quantitative Criminology, 18(1), 23-43.
#' @keywords aoristic
NULL
#' @name arlington
#' @title arlington burglary incident data (data frame)
#' @description A data frame of burglary incidents in the Arlington PD with four fields (DateTimeFrom, DateTimeTo, lon/lat).
#' @docType data
#' @usage data(aoristic)
#' @source Arlington PD 
#' @author George Kikuchi, 2013-09-13
NULL
#' @name CouncilDistrict
#' @title Council district (spatial polygon data frame)
#' @description A spatial polygon data frame of council districts in the Arlington PD jurisdiction.
#' @docType data
#' @usage data(aoristic)
#' @source City of Arlington 
#' @author George Kikuchi, 2013-09-13
NULL
#' @name aoristic
#' @title Sample data of crime (df) and council district (spdf)
#' @description Sample data of burglary incident data frame and spatial polygon data frame of council districts.
#' @docType data
#' @usage data(aoristic)
#' @source City of Arlington 
#' @author George Kikuchi, 2013-09-13
NULL