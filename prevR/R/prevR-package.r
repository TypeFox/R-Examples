#' Estimating regional trends of a prevalence from a DHS.
#'
#' \pkg{prevR} allows spatial estimation of a prevalence surface or a relative risks surface, 
#' using data from a Demographic and Health Survey (DHS) or an analog survey.
#'
#' @details \tabular{ll}{
#' Package: \tab prevR\cr
#' Type: \tab Package\cr
#' Licence: \tab CeCILL-C - \url{http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html}\cr
#' Website: \tab \url{http://www.ceped.org/prevR},\cr
#' \tab \url{http://joseph.larmarange.net/prevR}\cr
#' }
#' 
#' This package performs a methodological approach for spatial estimation of regional trends 
#' of a prevalence using data from surveys using a stratified two-stage sample design 
#' (as Demographic and Health Surveys). In these kind of surveys, positive and control cases
#' are spatially positioned at the centre of their corresponding surveyed cluster.
#' 
#' This package provides functions to estimate a prevalence surface using a kernel estimator 
#' with adaptative bandwiths of equal number of persons surveyed (a variant of the nearest 
#' neighbour technique) or with fixed bandwiths. The prevalence surface could also be calculated 
#' using a spatial interpolation (kriging or inverse distance weighting) after a moving average
#' smoothing based on circles of equal number of observed persons or circles of equal radius.
#' 
#' With the kernel estimator approach, it's also possible to estimate a surface of relative risks.
#' 
#' For a quick demo, enter \code{quick.prevR(fdhs)}.
#' 
#' For a full demo, enter \code{demo(prevR)}.
#' 
#' The content of \pkg{prevR} can be broken up as follows:
#' 
#' \emph{Datasets}\cr
#' \code{\link{fdhs}} is a fictive dataset used for testing the package.\cr
#' \code{\link{TMWorldBorders}} provides national borders of every countries in the World and 
#'   could be used to define the limits of the studied area.
#' 
#' \emph{Creating objects}\cr
#' \pkg{prevR} functions takes as input ojects of class \code{\link[=prevR-class]{prevR}}.\cr
#' \code{\link{import.dhs}} allows to import easily, through a step by step procedure, 
#'   data from a DHS (Demographic and Health Surveys) downloaded from 
#'   \url{http://www.measuredhs.com}.\cr
#' \code{\link{as.prevR}} is a generic function to create an object of class 
#'   \code{\link[=prevR-class]{prevR}}.\cr
#' \code{\link{create.boundary}} could be used to select borders of a country and 
#'   transfer them to \code{\link{as.prevR}} in order to define the studied area.
#' 
#' \emph{Data visualisation}\cr
#' Methods \code{\link[=show,prevR-method]{show}}, \code{\link[=print,prevR-method]{print}} 
#'   and \code{\link[=summary,prevR-method]{summary}} display a summary of a object of class 
#'   \code{\link[=prevR-class]{prevR}}.\cr
#' The method \code{\link[=plot,prevR-method]{plot}} could be used on a object of class 
#'   \code{\link[=prevR-class]{prevR}} for visualising the studied area, spatial position 
#'   of clusters, number of observations or number of positive cases by cluster.
#' 
#' \emph{Data manipulation}\cr
#' The method \code{\link[=changeproj,prevR-method]{changeproj}} changes the projection 
#'   of the spatial coordinates.\cr
#' The method \code{\link[=as.data.frame.prevR]{as.data.frame}} converts an object of 
#'   class \code{\link[=prevR-class]{prevR}} into a data frame.\cr
#' The method \code{\link[=export,prevR-method]{export}} export data and/or the studied 
#'   area in a text file, a dbf file or a shapefile.
#' 
#' \emph{Data analysis}\cr
#' \code{\link[=rings,prevR-method]{rings}} calculates rings of equal number of 
#'   observations and/or equal radius.\cr
#' \code{\link[=kde,prevR-method]{kde}} calculates a prevalence surface or a relative 
#'   risks surface using gaussian kernel density estimators (kde) with adaptative bandwiths.\cr
#' \code{\link[=krige,prevR-method]{krige}} executes a spatial interpolation using an 
#'   ordinary kriging.\cr
#' \code{\link[=idw,prevR-method]{idw}} executes a spatial interpolation using an inverse 
#'   distance weighting (idw) technique.
#' 
#' \emph{Results visualisation and export}\cr
#' Outputs of \code{\link[=kde,prevR-method]{kde}}, \code{\link[=krige,prevR-method]{krige}} 
#'   and \code{\link[=idw,prevR-method]{idw}} are objects of class 
#'   \code{\link[sp:SpatialPixelsDataFrame-class]{SpatialPixelsDataFrame}}\{\pkg{sp}\}.\cr
#' Results could be plotted using the function  \code{\link[sp]{spplot}}\{\pkg{sp}\}.\cr
#' \pkg{prevR} provides several continuous color palettes (see \code{\link{prevR.colors}}) 
#'   compatible with \code{\link[sp]{spplot}}.\cr
#' Calculated surfaces could be export using the function 
#'   \code{\link[maptools]{writeAsciiGrid}}\{\pkg{maptools}\}.
#' 
#' @author
#' Joseph Larmarange \email{joseph.larmarange@@ird.fr}\cr
#' IRD - CEPED (UMR 196 Université Paris Descartes Ined IRD)
#' 
#' @section Acknowledgement: \pkg{prevR} has been developed with funding from the French National Agency
#' for Research on AIDS and Viral Hepatitis (ANRS - \url{http://www.anrs.fr}) and the French Research
#' Institute for Development (IRD - \url{http://www.ird.fr}), and technical support from  LYSIS 
#' (info@@lysis-consultants.fr).
#' 
#' @section Citation: To cite \pkg{prevR}:\cr
#' Larmarange Joseph, Vallo Roselyne, Yaro Seydou, Msellati Philippe and Meda Nicolas (2011) "Methods for
#' mapping regional trends of HIV prevalence from Demographic and Health Surveys (DHS)", \emph{Cybergeo:
#' European Journal of Geography}, no 558, \url{http://cybergeo.revues.org/24606}, 
#' DOI: 10.4000/cybergeo.24606.
#' 
#' @references
#' Larmarange Joseph and Bendaud Victoria (2014) "HIV estimates at second subnational level 
#' from national population-based survey", \emph{AIDS}, n° 28, p. S469-S476, 
#' DOI: 10.1097/QAD.0000000000000480.
#' 
#' Larmarange Joseph, Vallo Roselyne, Yaro Seydou, Msellati Philippe and Meda Nicolas (2011) 
#' "Methods for mapping regional trends of HIV prevalence from Demographic and Health Surveys (DHS)",
#' \emph{Cybergeo: European Journal of Geography}, n° 558, \url{http://cybergeo.revues.org/24606}, 
#' DOI: 10.4000/cybergeo.24606.
#' 
#' Larmarange Joseph (2007) \emph{Prévalences du VIH en Afrique : validité d'une mesure}, 
#' PhD thesis in demography, directed by Benoît Ferry, université Paris Descartes, 
#' \url{http://tel.archives-ouvertes.fr/tel-00320283}.
#' 
#' Larmarange Joseph, Vallo Roselyne, Yaro Seydou, Msellati Philippe Meda Nicolas and 
#' Ferry Benoît (2006), "Cartographier les données des enquêtes démographiques et de santé 
#' à partir des coordonnées des zones d'enquête", \emph{Chaire Quételet, 
#' 29 novembre au 1er décembre 2006}, Université Catholique de Louvain, Louvain-la-Neuve, 
#' Belgique, \url{http://www.uclouvain.be/13881.html}.
#' @examples 
#' \dontrun{
#' par(ask = TRUE)
#' # Creating an object of class prevR
#' col <- c(id = "cluster", 
#'         x = "x",
#'         y = "y",
#'         n = "n",
#'         pos = "pos",
#'         c.type = "residence",
#'         wn = "weighted.n",
#'         wpos = "weighted.pos"
#' )
#' dhs <- as.prevR(fdhs.clusters,col, fdhs.boundary)
#' 
#' str(dhs)
#' print(dhs)
#' 
#' plot(dhs, main="Clusters position")
#' plot(dhs, type="c.type", main="Clusters by residence")
#' plot(dhs, type="count", main="Observations by cluster")
#' plot(dhs, type="flower", main="Positive cases by cluster")
#' 
#' # Changing coordinates projection
#' plot(dhs,axes=TRUE)
#' dhs <- changeproj(dhs,
#'                  "+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#' print(dhs)
#' plot(dhs, axes=TRUE)
#' 
#' # Calculating rings of equal number of observations for different values of N
#' dhs <- rings(dhs,N=c(100,200,300,400,500))
#' print(dhs)
#' summary(dhs)
#' 
#' # Prevalence surface for N=300
#' prev.N300 <- kde(dhs, N=300, nb.cells=200)
#' spplot(prev.N300, 'k.wprev.N300.RInf',
#'       cuts=100, col.regions=prevR.colors.red(101),
#'       main="Regional trends of prevalence (N=300)"
#' )
#' 
#' # Smoothing ring radii surface (spatial interpolation by kriging)
#' radius.N300 <- krige('r.radius', dhs, N=300, nb.cells=200)
#' spplot(radius.N300, 
#'       cuts=100, col.regions=prevR.colors.blue(101),
#'       main="Radius of circle (N=300)"
#' )
#' par(ask = FALSE)
#' }
#' @docType package
#' @name prevR-package
#' @keywords package
#' @import sp
#' @import rgdal
#' @import ggplot2
#' @import directlabels
#' @importFrom gstat idw krige vgm as.vgm.variomodel fit.variogram variogram
#' @importFrom fields rdist rdist.earth
#' @importFrom GenKern KernSur
#' @importFrom methods setClass setGeneric setMethod
#' @importFrom maptools writePointsShape writePolyShape
#' @importFrom foreign write.dbf read.dbf read.spss
#' @importFrom grDevices col2rgb dev.new gray hsv
#' @importFrom graphics legend points rect sunflowerplot text title
#' @importFrom methods as is new slot slot<- slotNames
#' @importFrom stats na.omit optim quantile
#' @importFrom utils alarm de edit menu select.list setTxtProgressBar txtProgressBar write.table
NULL

#' Fictitious data generated by a DHS simulation.
#'
#' Data set generated by a Demographic and Health Survey (DHS) simulation on a fictitious 
#' country with a national prevalence of 10\%, 8000 having been surveyed, distributed in 401 clusters. 
#' This dataset is composed of 3 objects:\itemize{
#'   \item \code{fdhs.clusters}: data frame (one line per cluster).
#'   \item \code{fdhs.boundary}: object of class 
#'     \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}} corresponding to the borders of the fictitious country.
#'   \item \code{fdhs}: object of class \code{\link[=prevR-class]{prevR}} 
#'     returned by \code{\link{as.prevR}} using the two previous objects.
#' }
#'
#' @examples 
#' \dontrun{
#'   str(fdhs)
#'   str(fdhs.clusters)
#'   str(fdhs.boundary)
#'   demo(prevR)
#' } 
#' @name fdhs
#' @aliases fdhs.boundary fdhs.clusters
#' @docType data
#' @keywords datasets
NULL

#' Dataset "TM World Borders Dataset 0.3".
#' 
#' This dataset provides boundaries of all countries in the world, in decimal degrees. 
#' Available variables are:\itemize{
#'   \item "FIPS" FIPS 10-4 Country Code.
#'   \item "ISO2" ISO 3166-1 Alpha-2 Country Code.
#'   \item "ISO3" ISO 3166-1 Alpha-3 Country Code.
#'   \item "UN" ISO 3166-1 Numeric-3 Country Code.
#'   \item "NAME" Name of country/area.
#'   \item "AREA" Land area, FAO Statistics (2002).
#'   \item "POP2005" Population, World Population Prospects (2005).
#'   \item "REGION" Macro geographical (continental region), UN Statistics.
#'   \item "SUBREGION" Geographical sub-region, UN Statistics.
#'   \item "LON" Longitude.
#'   \item "LAT" Latitude.
#' }
#' @format Object of class \code{\link[sp:SpatialPolygonsDataFrame-class]{SpatialPolygonsDataFrame}}.
#' @source Provided by Bjorn Sandvik on \url{http://thematicmapping.org/downloads/world_borders.php}.
#' The dataset was derived by Schuyler Erle from public domain sources. Sean Gilles did some clean up 
#' and made some enhancements. The dataset is available under a 
#' \emph{Creative Commons Attribution-Share Alike License} 
#' (\url{http://creativecommons.org/licenses/by-sa/3.0/}).
#' @note The boundaries, names designations used do not imply official endorsement or acceptance 
#' by the authors. Use this dataset with care, as several of the borders are disputed
#' @examples 
#' plot(TMWorldBorders)
#' @name TMWorldBorders
#' @docType data
#' @keywords datasets spatial
NULL