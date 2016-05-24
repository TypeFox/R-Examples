#' Atomic Force Microscopy images tools
#' 
#' The AFM package provides statistics analysis tools for Atomic Force Microscopy image analysis.\cr
#' Licence: Affero GPL v3
#' 
#' A graphical user interface is available by using \code{\link{runAFMApp}} command.
#' 
#' Several high level functions are :
#' \itemize{
#'   \item create your AFM image from a list of measured  heights (see example section of \code{\link{AFMImage}})
#'   \item import your image from Nanoscope Analysis (TM) tool (\code{\link{importFromNanoscope}}) 
#'   \item check if your sample is normally distributed and isotropic and get a pdf report (\code{\link{generateCheckReport}})
#'   \item perform variance (variogram), roughness against lengthscale, fractal analysis and get a pdf report (\code{\link{generateReport}})
#' }
#' 
#' Other functions are :
#' \itemize{
#'   \item check sample: for normality (\code{\link{checkNormality}}) and for isotropy (\code{\link{checkIsotropy}})
#'   \item calculate total RMS roughness: quick calculation of total root mean square roughness(\code{\link{totalRMSRoughness}})
#'   \item calculate omnidirectional variogram: calculate estimated variogram (\code{\link{calculateOmnidirectionalVariogram}})
#'   \item calculate roughness against lenghscale and Power Spectrum Density (PSD): calculate roughness against length scale (\code{\link{RoughnessByLengthScale}}), PSD 1D (\code{\link{PSD1DAgainstFrequency}}) or PSD 2D (\code{\link{PSD2DAgainstFrequency}}) against frequencies 
#'   \item calculate fractal dimension and scale: use (\code{\link{getFractalDimensions}}) function
#'   \item print in 3D (3D print) (\code{\link{exportToSTL}}) your AFM image
#' }
#' 
#' An EC2 instance is available for basic testing at the following address: \url{http://www.afmist.org}
#' 
#' Note: To use with a Brucker(TM) Atomic Force Microscope, use nanoscope analysis(TM) software and
#' \itemize{
#'   \item Use the "Flatten" function.
#'   \item Save the flattened image.
#'   \item Use the "Browse Data Files" windows, right click on image name and then Export the AFM image with the headers and the "Export> ASCII" contextual menu option. 
#' }
#' 
#' @examples
#' \dontrun{
#'   library(AFM)
#' # Analyse the AFMImageOfRegularPeaks AFM Image from this package
#'   data("AFMImageOfRegularPeaks")
#'   AFMImage<-AFMImageOfRegularPeaks
#' # exportDirectory="C:/Users/my_windows_login" or exportDirectory="/home/ubuntu"
#'   exportDirectory=tempdir()
#'   AFMImage@@fullfilename<-paste(exportDirectory,"AFMImageOfRegularPeaks.txt",sep="/")
#'   
#' # Start to check if your sample is normaly distributed and isotropic.
#'   generateCheckReport(AFMImage)
#'   
#' # If the sample is normaly distributed and isotropic, generate a full report
#'   generateReport(AFMImage)
#'   }
#' @references 
#' Gneiting2012, Tilmann Gneiting, Hana Sevcikova and Donald B. Percival 'Estimators of Fractal Dimension: Assessing the Roughness of Time Series and Spatial Data - Statistics in statistical Science, 2012, Vol. 27, No. 2, 247-277' \cr\cr
#' Olea2006, Ricardo A. Olea "A six-step practical approach to semivariogram modeling", 2006, "Stochastic Environmental Research and Risk Assessment, Volume 20, Issue 5 , pp 307-318" \cr\cr
#' Sidick2009, Erkin Sidick "Power Spectral Density Specification and Analysis of Large Optical Surfaces", 2009, "Modeling Aspects in Optical Metrology II, Proc. of SPIE Vol. 7390 73900L-1"
#' @seealso \code{\link{gstat}}, \code{\link{fractaldim}}, \code{\link{rgl}}
#' @author M.Beauvais, J.Landoulsi, I.Liascukiene
#' @docType package
#' @name AFM
#' @import data.table
#' @import gstat
#' @import rgl
#' @import ggplot2
#' @import png
#' @import plyr
#' @import methods
#' @import igraph
#' @import methods
#' @importFrom sp coordinates coordinates<- CRS is.projected proj4string proj4string<- spplot
#' @importFrom fftwtools fftw2d
#' @importFrom fractaldim fd.estim.filter1 fd.estim.isotropic fd.estim.squareincr fd.estim.transect.incr1 fd.estim.transect.var
#' @importFrom grid grid.layout grid.newpage grid.text pushViewport viewport gpar grid.raster
#' @importFrom gridExtra tableGrob ttheme_default
#' @importFrom grDevices blues9 dev.off heat.colors pdf png
#' @importFrom moments skewness
#' @importFrom pracma ceil meshgrid
#' @importFrom rglwidget renderRglwidget rglwidget sceneChange registerSceneChange rglwidgetOutput
#' @importFrom shiny actionButton downloadButton downloadHandler fileInput h3 hr htmlOutput HTML imageOutput isolate mainPanel navbarMenu navbarPage observeEvent plotOutput radioButtons reactive reactiveValues renderImage renderPlot renderTable renderUI shinyServer shinyUI sidebarLayout sidebarPanel sliderInput tableOutput tabPanel updateSliderInput uiOutput
#' @importFrom shinyjs disable enable useShinyjs
#' @importFrom stats coefficients cor dist dnorm lm na.omit sd var
#' @importFrom stringr str_sub str_replace_all
#' @importFrom utils head installed.packages read.table tail write.table
NULL