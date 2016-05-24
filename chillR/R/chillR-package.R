

#' Weather stations in California
#' 
#' This is a list of weather stations in California that are contained in the
#' UC IPM database. This can also be generated with
#' make_california_UCIPM_station_list(), but this takes quite a while. So this
#' dataset is supposed to be a shortcut to this.
#' 
#' 
#' @name california_stations
#' @docType data
#' @format a data.frame containing stations from the California UC IPM database
#' (), with the following columns: "Name", "Code", "Interval", "Lat", "Long",
#' "Elev".  \describe{ \item{list("Name")}{mame of the weather station}
#' \item{list("Code")}{ code of the weather station, indicating the name and
#' the database it comes from} \item{list("Interval")}{ period of available
#' data (as character string)} \item{list("Lat")}{ latitude of the station}
#' \item{list("Long")}{ longitude of the station} \item{list("Elev")}{
#' elevation of the station} }
#' @source UC IPM website: http://www.ipm.ucdavis.edu/WEATHER/index.html
#' @keywords datasets
#' @importFrom pls plsr explvar crossval
#' @importFrom Kendall Kendall
#' @importFrom fields Krig surface tim.colors image.plot predictSurface
#' @importFrom sp spDistsN1
#' @importFrom readxl read_excel
#' @importFrom XML htmlParse getNodeSet xmlToDataFrame getChildrenStrings
#' @importFrom httr POST content
#' @importFrom grDevices bmp colorRampPalette dev.off png tiff gray.colors rainbow rgb
#' @importFrom graphics arrows axis box grconvertX lines mtext par plot points barplot contour image layout rect text
#' @importFrom stats aggregate coef sd median
#' @importFrom utils stack write.csv download.file read.csv read.fwf read.table unzip

#' @examples
#' 
#' data(california_stations)
#' 
NULL





#' Statistical methods for phenology analysis in temperate fruit trees
#' 
#' The phenology of plants (i.e. the timing of their annual life phases)
#' depends on climatic cues. For temperate trees and many other plants, spring
#' phases, such as leaf emergence and flowering, have been found to result from
#' the effects of both cool (chilling) conditions and heat. Fruit tree
#' scientists (pomologists) have developed some metrics to quantify chilling
#' and heat. 'chillR' contains functions for processing temperature records
#' into chilling (Chilling Hours, Utah Chill Units and Chill Portions) and heat
#' units (Growing Degree Hours). Regarding chilling metrics, Chill Portions are
#' often considered the most promising, but they are difficult to calculate.
#' This package makes it easy. 'chillR' also contains procedures for conducting
#' a PLS analysis relating phenological dates (e.g. bloom dates) to either mean
#' temperatures or mean chill and heat accumulation rates, based on long-term
#' weather and phenology records.
#' 
#' \tabular{ll}{ Package: \tab chillR\cr Type: \tab Package\cr License: \tab The "GNU General Public
#' License" version 3\cr }
#' 
#' @name chillR-package
#' @aliases chillR-package chillR
#' @docType package
#' @author Eike Luedeling Maintainer: Eike Luedeling <eike@@eikeluedeling.com>
#' @references Applications of some of the methods in the package:
#' 
#' Luedeling E, Zhang M, Luedeling V and Girvetz EH, 2009. Sensitivity of
#' winter chill models for fruit and nut trees to climatic changes expected in
#' California's Central Valley. Agriculture, Ecosystems and Environment 133,
#' 23-31
#' 
#' Luedeling E, Zhang M, McGranahan G and Leslie C, 2009. Validation of winter
#' chill models using historic records of walnut phenology. Agricultural and
#' Forest Meteorology 149, 1854-1864
#' 
#' Luedeling E and Brown PH, 2011. A global analysis of the comparability of
#' winter chill models for fruit and nut trees. International Journal of
#' Biometeorology 55, 411-421
#' 
#' Luedeling E, Kunz A and Blanke M, 2011. Mehr Chilling fuer Obstbaeume in
#' waermeren Wintern? (More winter chill for fruit trees in warmer winters?).
#' Erwerbs-Obstbau 53, 145-155
#' 
#' Luedeling E, Guo L, Dai J, Leslie C, Blanke M, 2013. Differential responses
#' of trees to temperature variation during the chilling and forcing phases.
#' Agricultural and Forest Meteorology 181, 33-42.
#' 
#' Review on chilling models in a climate change context:
#' 
#' Luedeling E, 2012. Climate change impacts on winter chill for temperate
#' fruit and nut production: a review. Scientia Horticulturae 144, 218-229
#' 
#' The chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords package
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2004&!(
#'    KA_weather$Year==2005&KA_weather$Month<6)),])
#' 
#' PLS_results<-PLS_pheno(
#'   weather_data=KA_weather,
#'   split_month=6,   #last month in same year
#'   bio_data=KA_bloom)
#'   
#' PLS_results_path<-paste(getwd(),"/PLS_output",sep="")
#'   
#' plot_PLS(PLS_results,PLS_results_path)
#' 
#' stack<-stack_hourly_temps(weather,latitude=50.4)
#' cc<-chilling(stack,305,60)
#' 
NULL





#' Cherry bloom data for Klein-Altendorf, Germany
#' 
#' Bloom data of sweet cherry var. 'Schneiders spaete Knorpelkirsche' recorded
#' at Klein-Altendorf, Germany, the experimental station of the University of
#' Bonn
#' 
#' 
#' @name KA_bloom
#' @docType data
#' @format A data frame with the following 2 variables.  \describe{
#' \item{Year}{a numeric vector, indicating the observation year}
#' \item{pheno}{ a vector that, when coerced by as.numeric, contains
#' bloom data in Julian dates (day of the year)} }
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' @source data were collected by Achim Kunz and Michael Blanke, University of
#' Bonn
#' @keywords datasets
#' @examples
#' 
#' data(KA_bloom)
#' 
NULL





#' Weather data for Klein-Altendorf, Germany
#' 
#' Daily temperature data from Klein-Altendorf, Germany, for use in combination
#' with the example phenology dataset KA_bloom.
#' 
#' 
#' @name KA_weather
#' @docType data
#' @format A data frame with observations on the following 5 variables.
#' \describe{ \item{Year}{a numeric vector - the observation year}
#' \item{Month}{a numeric vector - the observation month}
#' \item{Day}{a numeric vector - the observation day}
#' \item{Tmax}{a numeric vector - daily maximum temperature}
#' \item{Tmin}{a numeric vector - daily minimum temperature} }
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' @source data were collected by Achim Kunz and Michael Blanke, University of
#' Bonn
#' @keywords datasets
#' @examples
#' 
#' data(KA_weather)
#' 
NULL



