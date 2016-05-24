#https://github.com/variani/pckdev/wiki/Documenting-with-roxygen2
#http://stackoverflow.com/questions/7356120/how-to-properly-document-s4-methods-using-roxygen2

#' Histogram-Valued Data Analysis
#' 
#' We consider histogram-valued data, i.e., data described by univariate histograms. The
#' methods and the basic statistics for histogram-valued data are mainly based
#' on the L2 Wasserstein metric between distributions, i.e., a Euclidean metric
#' between quantile functions. The package contains unsupervised classification techniques, 
#' least square regression and tools for histrogram-valued data and for histogram time series.
#'
#' 
#' \tabular{ll}{ Package: \tab HistDAWass\cr Type: \tab Package\cr Version:
#' \tab 0.1.1\cr Date: \tab 2014-09-17\cr License: \tab GPL (>=2)\cr Depends:
#' \tab methods\cr } ~~ An overview of how to use the package, including the
#' most important functions ~~
#' 
#' @name HistDAWass-package
#' @aliases HistDAWass-package HistDAWass
#' @docType package
#' @author Antonio Irpino <antonio.irpino@@unina2.it>
#' @references
#' Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI 10.1007/s11634-014-0176-4\cr
#' #' Irpino, A. and Romano, E. (2007): \emph{Optimal histogram representation of large data sets: 
#' Fisher vs piecewise linear approximations}. RNTI E-9, 99-110.
#' @import ggplot2 grid
#' @keywords package
#' @examples
#' 
#' # Generating a list of distributions
#' a<-vector("list",4)
#' a[[1]]<-distributionH(x=c(80,100,120,135,150,165,180,200,240),
#'                       p=c(0,0.025,0.1,0.275,0.525,0.725,0.887,0.975,1))
#' a[[2]]<-distributionH(x=c(80,100,120,135,150,165,180,195,210,240),
#'                       p=c(0,0.013,0.101,0.255,0.508,0.718,0.895,0.961,0.987,1))
#' a[[3]]<-distributionH(x=c(95,110,125,140,155,170,185,200,215,230,245),
#'                       p=c(0,0.012,0.041,0.154,0.36,0.595,0.781,0.929,0.972,0.992,1))
#' a[[4]]<-distributionH(x=c(105,120,135,150,165,180,195,210,225,240,260),
#'                       p=c(0,0.009,0.035,0.081,0.186,0.385,0.633,0.832,0.932,0.977,1))
#' # Generating a list of names of observations
#' namerows<-list( 'u1'  , 'u2') 
#' # Generating a list of names of variables
#' namevars<-list( 'Var_1'  , 'Var_2') 
#' # creating the MatH  
#' Mat_of_distributions<-MatH(x=a, nrows = 2, ncols = 2, 
#'                             rownames=namerows, varnames=namevars, by.row=FALSE )
#' 
NULL


#' Full Ozone dataset for Histogram data analysis
#' 
#' The dataset contains MatH (matrix of histogram-valued data) object This data
#' set list 78 stations located in the USA recording four variables, without
#' missing data.
#' 
#' 
#' @name OzoneFull
#' @docType data
#' @format a \code{MatH} istance, 1 row per station.
#' @author Antonio Irpino, 2014-10-05
#' @source http://java.epa.gov/castnet/epa_jsp/prepackageddata.jsp
#' ftp://ftp.epa.gov/castnet/data/metdata.zip
NULL



#' Complete Ozone dataset for Histogram data analysis
#' 
#' The dataset contains MatH (matrix of histogram-valued data) object This data
#' set list 84 stations located in the USA recording four variables. Some
#' stations contains missing data.
#' 
#' 
#' @name OzoneH
#' @docType data
#' @format a \code{MatH} istance, 1 row per station.
#' @author Antonio Irpino, 2014-10-05
#' @source http://java.epa.gov/castnet/epa_jsp/prepackageddata.jsp
#' ftp://ftp.epa.gov/castnet/data/metdata.zip
NULL

#' Blood dataset for Histogram data analysis
#' 
#' The dataset contains a MatH (matrix of histogram-valued data) object This
#' data set list 14 groups of patients described by 3 variables.
#' 
#' 
#' @name BLOOD
#' @docType data
#' @format a \code{MatH} istance, 1 row per group.
#' @author Antonio Irpino, 2014-10-05
#' @source Billard L. and Diday E. (2006). Symbolic Data Analysis: Conceptual
#' Statistics and Data Mining, Wiley.
NULL

#' Blood dataset from Brito P. for Histogram data analysis
#' 
#' The dataset contains a MatH (matrix of histogram-valued data) object This
#' data set list 10 patients described by 2 variables.
#' 
#' 
#' @name BloodBRITO
#' @docType data
#' @format a \code{MatH} istance, 1 row per patient.
#' @author Antonio Irpino, 2014-10-05
#' @source Dias, S. and Brito P. Distribution and Symmetric Distribution
#' Regression Model for Histogram-Valued Variables, ArXiv, arXiv:1303.6199
#' [stat.ME]
NULL

#' Age pyramids of all the countries of the World in 2014
#' 
#' The dataset contains a MatH (matrix of histogram-valued data) object, with three 
#' hisogram-valued variables, the 5-years age (relative frequencies) distribution of all the population, 
#' of the male and of the female population of 228 countries of the World. The  first row is the World data.
#' Thus it contains 229 rows(228 countries plus the World) and 3 variables: "Both.Sexes.Population", 
#' "Male.Population", "Female.Population" 
#' 
#' 
#' @name Age_Pyramids_2014
#' @docType data
#' @format a \code{MatH} object, a matrix of distributions.
#' @author Antonio Irpino, 2014-10-05
#' @source United States Census Bureau \url{http://www.census.gov/data.html}
NULL

#' Agronomique data
#' 
#' A dataset with the distributions of marginal costs of farms in 22 France regions. It contains four 
#' histogram variables: "Y_TSC" (Total costs of a farm),  "X_Wheat" (Costs for Wheat), "X_Pig" (Costs for Pigs)
#' "X_Cmilk" (Costs for Cow Milk)
#' 
#' @name Agronomique
#' @docType data
#' @format a \code{MatH} object, a matrix of distributions.
#' @author Antonio Irpino, 2014-10-05
#' @source Rosanna Verde, Antonio Irpino, Second University of Naples; Dominique Desbois, UMR Economie
#' publique, INRA-AgroParisTech, How to cope with modelling and privacy concerns? A regression model and a visualization tool for
#' aggregated data,  Conference of European Statistics Stakeholders,
#' Rome, November, 24-25,2014
NULL

#' A monthly climatic dataset of China
#' 
#' A dataset with the distributions of some climatic variables collected for each month in 60 stations of China.
#' The collected variables are 168 i.e. 14 climatic variables observed for 12 months. The 14 variables are the following:
#'  mean station pressure (mb), mean temperature, mean maximum temperature, mean minimum temperature,
#'   total precipitation (mm), sunshine duration (h), mean cloud amount (percentage of sky cover),
#'   mean relative humidity (%), snow days (days with snow cover), dominant wind direction (degrees),
#'    mean wind speed (m/s), dominant wind frequency (%), extreme maximum temperature (C),
#'    extreme minimum temperature. 
#' Use the command \code{ get.MatH.main.info(China_Month)} for rapid info.
#'
#' @name China_Month
#' @docType data
#' @format a \code{MatH} object, a matrix of distributions.
#' @author Antonio Irpino, 2014-10-05
#' @source raw data are available here: \url{http://cdiac.ornl.gov/ndps/tr055.html}
NULL

#' A seasonal climatic dataset of China
#' 
#' A dataset with the distributions of some climatic variables collected for each season in 60 stations of China.
#' The collected variables are 56 i.e. 14 climatic variables observed for 4 seasons. The 14 variables are the following:
#'  mean station pressure (mb), mean temperature, mean maximum temperature, mean minimum temperature,
#'   total precipitation (mm), sunshine duration (h), mean cloud amount (percentage of sky cover),
#'   mean relative humidity (%), snow days (days with snow cover), dominant wind direction (degrees),
#'    mean wind speed (m/s), dominant wind frequency (%), extreme maximum temperature (C),
#'    extreme minimum temperature. 
#' Use the command \code{ get.MatH.main.info(China_Seas)} for rapid info.
#'
#' @name China_Seas
#' @docType data
#' @format a \code{MatH} object, a matrix of distributions.
#' @author Antonio Irpino, 2014-10-05
#' @source raw data are available here: \url{http://cdiac.ornl.gov/ndps/tr055.html}. 
#' Climate Data Bases of the People's Republic of China 1841-1988 (TR055)
#' DOI: 10.3334/CDIAC/cli.tr055
NULL

#' Stations coordinates of China_Month and China_Seas datasets
#' 
#' A dataset containing the geographical coordinates of stations described in China_Month and China_Seas datasets
#'
#' 
#' @name stations_coordinates
#' @docType data
#' @format a data.frame
#' @author Antonio Irpino, 2014-10-05
#' @source raw data are available here: \url{http://cdiac.ornl.gov/ndps/tr055.html}. 
#' Climate Data Bases of the People's Republic of China 1841-1988 (TR055)
#' DOI: 10.3334/CDIAC/cli.tr055
NULL

#' A histogram-valued dataset of returns
#' 
#' A histogram-valued dataset of returns of dollar vs yen change rates
#' 
#' @name RetHTS
#' @docType data
#' @format a \code{MatH} object, a matrix of distributions.
#' @author Antonio Irpino, 2014-10-05
NULL















 






 




