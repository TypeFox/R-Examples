#' Distance sampling
#'
#' \code{Distance} is a simple way to fit detection functions and estimate abundance using distance sampling methodology.
#'
#' Underlying \code{Distance} is the package \code{mrds}, for more advanced analyses (such as those involving double observer surveys) one may find it necessary to use \code{mrds}.
#'
#' Further information on distance sampling methods and example code is available at \url{http://distancesampling.org/R/}.
#'
#' For help with distance sampling and this package, there is a Google Group \url{https://groups.google.com/forum/#!forum/distance-sampling}.
#'
#' @name Distance-package
#' @import mrds
#' @aliases Distance-package Distance
#' @docType package
#' @author David L. Miller <dave@@ninepointeightone.net>
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#'
#' Marques, F.F.C. and S.T. Buckland. 2004. Covariate models for the detection
#'   function. In: Advanced Distance Sampling, eds. S.T. Buckland,
#'   D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas.
#'   Oxford University Press.
#'
#' @keywords Statistical Models
#'
NULL


#' Simulated minke whale data
#'
#' Data simulated from models fitted to 1992/1993 Southern Hemisphere minke whale data collected by the International Whaling Commission. See Branch and Butterworth (2001) for survey details (survey design is shown in figure 1(e)). Data simulated by David Borchers.
#'
#' Data are included here as both R data and as an Excel spreadsheet to illustrate the "flat file" input method. See \code{\link{flatfile}} for how to load this data and an example analysis.
#'
#' @references Branch, T.A. and D.S. Butterworth (2001) Southern Hemisphere minke whales: standardised abundance estimates from the 1978/79 to 1997/98 IDCR-SOWER surveys. Journal of Cetacean Research and Management 3(2): 143-174
#'
#' Hedley, S.L., and S.T. Buckland. Spatial Models for Line Transect Sampling. Journal of Agricultural, Biological, and Environmental Statistics 9, no. 2 (2004): 181-199. doi:10.1198/1085711043578.
#'
#' @name minke
#' @keywords datasets
#' @source Shipped with the DISTANCE Windows application.
#' @docType data
#' @format \code{data.frame} with 99 observations of 5 variables:
#' \tabular{ll}{\code{Region.Label} \tab stratum label (\code{"North"} or \code{"South"}) \cr
#' \code{Area} \tab stratum area\cr
#' \code{Sample.Label} \tab transect identifier\cr
#' \code{Effort} \tab transect length\cr
#' \code{distance} \tab observed distance}
#' @examples
#' data(minke)
#' head(minke)
NULL

#' Flat files
#'
#' \code{Distance} allows loading data as a "flat file" and analyse data (and obtain abundance estimates) straight away, provided that the format of the flat file is correct. One can provide the file as, for example, an Excel spreadsheet using \code{read.xls} in \pkg{gdata} or CSV using \code{\link{read.csv}}.
#'
#' Each row of the data table corresponds to one observation and must have a the following columns:
#' \tabular{ll}{\code{distance} \tab observed distance to object \cr
#' \code{Sample.Label} \tab Identifier for the sample (transect id) \cr
#' \code{Effort} \tab effort for this transect (e.g. line transect length or number of times point transect was visited) \cr
#' \code{Region.Label} \tab code for regional strata (see below) \cr
#' \code{Area} \tab area of the strata}
#'
#' Note that in the simplest case (one area surveyed only once) there is only one \code{Region.Label} and a single corresponding \code{Area} duplicated for each observation.
#'
#' The example given below was provided by Eric Rexstad.
#'
#' @name flatfile
#' @docType methods
#' @examples
#' \dontrun{
#' library(Distance)
#' # Need to have the gdata library installed from CRAN, requires a system
#' # with perl installed (usually fine for Linux/Mac)
#' library(gdata)
#'
#' # Need to get the file path first
#' # Going to the path given in the below, one can examine the format
#' minke.filepath <- system.file("minke.xlsx",package="Distance")
#'
#' # Load the Excel file, note that header=FALSE and we add column names after
#' minke <- read.xls(minke.filepath, stringsAsFactor=FALSE,header=FALSE)
#' names(minke) <- c("Region.Label", "Area", "Sample.Label", "Effort","distance")
#' # One may want to call edit(minke) or head(minke) at this point
#' # to examine the data format
#'
#' # Due to the way the file was saved and the default behaviour in R
#' # for numbers stored with many decimal places (they are read as strings
#' # rather than numbers, see str(minke)). We must coerce the Effort column
#' # to numeric
#' minke$Effort <- as.numeric(minke$Effort)
#'
#' ## perform an analysis using the exact distances
#' pooled.exact <- ds(minke, truncation=1.5, key="hr", order=0)
#' summary(pooled.exact)
#'
#'
#' ## Try a binned analysis
#' # first define the bins
#' dist.bins <- c(0,.214, .428,.643,.857,1.071,1.286,1.5)
#' pooled.binned <- ds(minke, truncation=1.5, cutpoints=dist.bins, key="hr", order=0)
#'
#' # binned with stratum as a covariate
#' minke$stratum <- ifelse(minke$Region.Label=="North", "N", "S")
#' strat.covar.binned <- ds(minke, truncation=1.5, key="hr",
#'                          formula=~as.factor(stratum), cutpoints=dist.bins)
#'
#' # Stratified by North/South
#' full.strat.binned.North <- ds(minke[minke$Region.Label=="North",],
#'                   truncation=1.5, key="hr", order=0, cutpoints=dist.bins)
#' full.strat.binned.South <- ds(minke[minke$Region.Label=="South",],
#'                      truncation=1.5, key="hr", order=0, cutpoints=dist.bins)
#'
#' ## model summaries
#' model.sel.bin <- data.frame(name=c("Pooled f(0)", "Stratum covariate",
#'                                    "Full stratification"),
#'                             aic=c(pooled.binned$ddf$criterion,
#'                                   strat.covar.binned$ddf$criterion,
#'                                   full.strat.binned.North$ddf$criterion+
#'                                   full.strat.binned.South$ddf$criterion))
#'
#' # Note model with stratum as covariate is most parsimonious
#' print(model.sel.bin)
#' }
NULL

