#' Descriptive statistics of the cumulative distribution function of a
#' continuous variable
#' 
#' This function returns summary statistics of the cumulative distribution
#' function of a continuous variable estimated with \pkg{spsurvey}-package.
#' 
#' The function \code{cont.analysis()} of \pkg{spsurvey}-package estimates the
#' population total, mean, variance, and standard deviation of a continuous
#' variable. It also estimates the standard error and confidence bounds of
#' these population estimates. In some cases it may be interesting to see all
#' estimates, for which one uses \code{all = TRUE}. However, in other
#' circumstances there might be interest only in taking a look at the estimated
#' population mean and standard deviation. Then the argument \code{all} has to
#' be set to \code{FALSE}.
#' 
#' @param obj Object containing the estimated cumulative distribution function
#' of the continuous variable. The resulting object of \code{cont.analysis()}
#' of \pkg{spsurvey}-package.
#' @param ind Indicator variable. The name of the continuous variable as
#' displayed in the resulting object of \code{cont.analysis()}.
#' @param all Summary statistics to be returned. The default option (\code{all
#' = TRUE}) returns all summary statistics available. If \code{all = FALSE},
#' then only estimated population mean and standard deviation are returned. See
#' \sQuote{Details}.
#' @return A \code{data.frame} containing summary statistics of the cumulative
#' distribution function of a continuous variable.
#' @author Alessandro Samuel-Rosa <\email{alessandrosamuelrosa@@gmail.com}>
#' @seealso \code{\link[spsurvey]{cont.analysis}}.
#' @references Kincaid, T. M. and Olsen, A. R. (2013). spsurvey: Spatial Survey
#' Design and Analysis. R package version 2.6. URL:
#' <http://www.epa.gov/nheerl/arm/>.
#' @keywords methods print
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Estimate the CDF
#' my.cdf <- spsurvey::cont.analysis(spsurvey.obj = my.spsurvey)
#' 
#' ## See indicator levels in the resulting object
#' levels(my.cdf$Pct$Indicator)
#' 
#' ## Return all summary statistics of indicator variable 'dx'
#' cdfStats(my.cdf, "dx", all = TRUE)
#' }
#' 
# FUNCTION #####################################################################
cdfStats <- 
  function(obj, ind, all = TRUE) {
    stats <- data.frame(obj$Pct[obj$Pct$Indicator == ind, 4:9][8:10, ],
                        row.names = NULL)
    if(all) {
      res <- stats  
    } else {
      res <- stats[1, 3]
    }
    res
  }
# End!
