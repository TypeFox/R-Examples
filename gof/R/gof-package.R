##' Model-diagnostics based on cumulative residuals
##' 
##' @name gof-package
##' @import utils stats
##' @aliases gof-package gof
##' @docType package
##' @author Klaus K. Holst <kkho@@biostat.ku.dk>
##'@seealso \code{\link[timereg]{cox.aalen}} in the \code{timereg}-package for
##' similar GoF-methods for survival-data.
##'@references D.Y. Lin and L.J. Wei and Z. Ying (2002) \emph{Model-Checking
##' Techniques Based on Cumulative Residuals}. Biometrics, Volume 58, pp 1-12.
##'
##' John Q. Su and L.J. Wei (1991) \emph{A lack-of-fit test for the mean function
##' in a generalized linear model}. Journal. Amer. Statist. Assoc., Volume 86,
##' Number 414, pp 420-426.
##' @useDynLib gof
##' @import utils
##' @keywords package
##' @examples
##'
##' example(cumres)
##'
NULL


##' Surgical Unit Data
##'
##' Surgical Unit Data used in the paper by Lin et al. (2002).  Survival time and
##' covariates for 54 patients undergoing liver surgery.
##'
##'
##' @name surgunit
##' @docType data
##' @references D.Y. Lin and L.J. Wei and Z. Ying (2002) \emph{Model-Checking
##' Techniques Based on Cumulative Residuals}. Biometrics, Volume 58, pp 1-12.
##' @source Neter, J., Kutner, M. H., Nachtsheim, C. J., and Wasserman, W.
##' (1996), Applied Linear Statistical Models, 4th edition. Chicago: Irwin
NULL



