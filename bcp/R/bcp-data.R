#' @title Array CGH data set of Coriell cell lines
#' @name coriell
#' @description These are two data array CGH studies of Coriell cell lines taken from the reference below.
#' @docType data
#' @usage coriell
#' @format A data frame containing five variables: first is clone name, second is clone chromosome, third is   clone position, fourth and fifth  are log2ratio for two cell lines.
#' @references 
#' \enumerate{ 
#' \item Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004), Circular binary segmentation for the analysis of array-based DNA copy number data, \emph{Biostatistics}, \bold{5}, 557-572. url: \url{http://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html}.
#' \item Snijders \emph{et al.} (2001), Assembly of microarrays for genome-wide measurement of DNA copy number, \emph{Nature Genetics}, \bold{29}, 263-264. 
#' }
#' @source \url{http://www.nature.com/ng/journal/v29/n3/full/ng754.html}
#' @keywords datasets
#' @examples 
#' demo(coriell)
#' 
NULL

#' @title Milling machine indentation data
#' @name lombard
#' @format A vector of length 100 containing the individual radii.
#' @docType data
#' @usage lombard
#' @description Radii of 100 circular indentations cut by a milling machine.
#' @references 
#' \enumerate{ 
#' \item Daniel Barry and J. A. Hartigan (1993), A Bayesian Analysis for Change Point Problems, \emph{Journal of The American Statistical Association}, \bold{88}, 309-19.
#' \item F. Lombard (1987), Rank Tests for Changepoint Problems, \emph{Biometrika}, \bold{74}, 615-624. 
#' }
#' @source The data is available online in the data archive of the
#' \emph{Journal of Applied Econometrics}. url: 
#'  \url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.
#' @keywords datasets
#' @examples  
#' \dontrun{
#' data(lombard)
#' # univariate change point analysis
#' bcp.model <- bcp(lombard, burnin=500, mcmc=5000, return.mcmc=TRUE)
#' 
#' # linear model change point analysis 
#' bcpr.model <- bcp(lombard, cbind(1:100), burnin=500, mcmc=5000, return.mcmc=TRUE)
#' 
#' plot(bcp.model, main="Lombard Milling Data")
#' plot(bcpr.model, main="Lombard Milling Data (with Regression Model)")
#' }
NULL


#' @title New Haven housing data
#' @name NewHavenHousing
#' @format A matrix containing location, 2011 assessed value, and physical characteristics in the columns.
#' @docType data
#' @usage NewHavenHousing
#' @description Location, 2011 assessed value, and physical characteristics for 244 houses in a region of New Haven, CT.
#' @source The data can be scraped from the New Haven, CT Online Assessment Database \url{http://gis.vgsi.com/newhavenct/}
#' @references 
#' Xiaofei Wang and John W. Emerson (2015). Bayesian Change Point Analysis of Linear Models on General Graphs, \emph{Working Paper}.
#' @keywords datasets
#' @examples  
#' \dontrun{
#' demo("NewHaven")
#' }
#' 
NULL


#' @title Quebec river streamflow data
#' @name QuebecRivers
#' @format A matrix containing streamflow amounts for the years 1972 to 1994.
#' @docType data
#' @usage QuebecRivers
#' @description Annual January to June streamflow amounts (measured in \eqn{L/(km^2 s)}) for four rivers in Quebec: Baleine, Churchill Falls, Manicouagan, and Romaine.
#' @source The data can be obtained from the Centre d'expertise hydrique Quebec. \url{https://www.cehq.gouv.qc.ca/hydrometrie/index-en.htm}
#' @references 
#' \enumerate{ 
#' \item Xiaofei Wang and John W. Emerson (2015). Bayesian Change Point Analysis of Linear Models on General Graphs, \emph{Working Paper}.
#' \item L. Perrault \emph{et al.} (2000). Retrospective multivariate Bayesian change-point analysis: a simultaneous single change in the mean of several hydrological sequences, \emph{Stochastic Environmental Research and Risk Assessment}, \bold{14}, 243-261.
#' }
#' @keywords datasets
#' @examples  
#' data("QuebecRivers")
#' bcpr.rivers <- bcp(QuebecRivers)
#' plot(bcpr.rivers, main="Quebec River Streamflow Change Point Analysis", 
#'      xlab="Year", xaxlab = 1972:1994)
#' 
NULL


#' @title US Ex-post Real Interest Rate data, 1961(1):1986(3)
#' @name RealInt
#' @format A quarterly time series from 1961(1) to 1986(3).
#' @docType data
#' @usage RealInt
#' @description US ex-post real interest rate: the three-month treasury bill
#' deflated by the CPI inflation rate.
#' @references 
#' \enumerate{ 
#' \item J. Bai and P. Perron (2003), Computation and Analysis of Multiple Structural Change Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22. \url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.
#' \item Achim Zeileis, Friedrich Leisch, Kurt Hornik, Christian Kleiber (2002), strucchange: An R Package for Testing for Structural Change in Linear Regression Models, \emph{Journal of Statistical Software}, \bold{7}(2), 1--38. \url{http://www.jstatsoft.org/v07/i02/}.
#' }
#' @source The data is available online in the data archive of the
#' \emph{Journal of Applied Econometrics}. url: 
#'  \url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.
#' @keywords datasets
#' @examples  
#' demo(RealInt)
#' 
NULL
