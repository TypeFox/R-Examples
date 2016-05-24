#' Optigrab  
#' 
#' @description 
#' The optigrab packages providesa function \code{opt_grab} to retrieve 
#' options/arguments from the command line.  It is useful for running R in 
#' batch mode with \code{R CMD BATCH ...} or \code{Rscript}.  
#' 
#' GNU-, Java- and Microsoft-style command line options are supported. GNU-style
#' is the default.
#' 
#' See the \emph{optigrab} vignettes or the github README file for details.
#'
#' 
#' @seealso 
#'   \code{\link{opt_get}} \cr
#'   \code{\link[base]{commandArgs}} \cr
#' 
#' @references 
#'   The Jerk. Dir. Carl Reiner. Perf Steve Martin, Bernadette Peters, 
#'   Caitlin Adams. Universal Pictures, 1979. \cr
#' 
#'   \url{http://www.gnu.org/prep/standards/standards.html} \cr
#' 
#'   \url{ https://github.com/gsscoder/commandline} \cr
#'   
#' @examples
#' 
#' \dontrun{ 
#'   opt_get( c("foo","f"))
#' }
#' 
#'   opts <- c( "--flag", "bar" ) 
#'   flag <- opt_get( c("foo","f"), opts=opts )  # bar
#'
#' @name optigrab
#' @docType package
#' @aliases optigrab-package 
#' @keywords package

NULL
