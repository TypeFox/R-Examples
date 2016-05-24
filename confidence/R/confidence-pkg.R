#'  Confidence Estimation of Environmental State Classifications
#' 
#'  This package can be used to estimate the confidence of state classifications (e.g., with the classification `bad', `moderate', `good') produced using environmental indicators and associated targets. The implementation closely follows Baggelaar et al. (2010) where the confidence intervals for the estimated multiyear averages are derived by assuming a Student's t distribution for the errors. For more information about the package see the package-vignette (type: \code{vignette("confidence")} at the \R-prompt to view the vignette.).
#'  
#'  @seealso \code{\link{conf}} and the package vignette mentioned above.
#'  
#'  @author Willem M.G.M. van Loon and Dennis J.J. Walvoort
#'    
#'  @references Baggelaar, P., O. van Tongeren, R. Knoben, and W. van Loon, 2010. Rapporteren van de betrouwbaarheid van KRW-beoordelingen (in Dutch, English translation: Reporting the confidence of WFD-assessments). H2O 16: 21--25
#' 
#'  @import knitr
#'  @import markdown
#'  @import plyr
#'  @import xtable
#'  @import ggplot2
#'  @import tcltk
#'  @name confidence-pkg
#'  @aliases confidence
#'  @docType package
#'  @keywords package
NULL