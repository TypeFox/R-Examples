#' Plot Statistical Power Curve
#' 
#' Generic \code{plot} function for the \code{pwrrasch} object, which 
#' plots the statistical power curve relating statistical power to sample size 
#'
#' @details Graphical parameters are:
#' \itemize{
#'  \item{\code{type}} The following values are possible: \code{"p"} for points,
#'                     \code{"l"} for lines, \code{"b"} for both point and lines
#'  \item{\code{pch}} see \link[graphics]{points}
#'  \item{\code{lty}} Line types can  be specified as an integer (\code{0} = blank, \code{1} = solid,
#'                    \code{2} = dashed, \code{3} = dotted, \code{4} = dotdash, \code{5} = longdash, 
#'                    \code{6} = twodash)
#'  \item{\code{lwd}} Positive numbers indicating line widths   
#'  \item{\code{legend}} Either the x and y coordinates to be used to position the legend or
#'                       keyword from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, 
#'                       \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} 
#'                       and \code{"center"}
#'  \item{\code{bty}} Allowed values are "o" (draw box around legend) and "n" (do not draw box around legend).                                                            
#' }
#' 
#' @param x                \code{pwrrasch} object.
#' @param plot.sig.level   If \code{TRUE}, nominal significance level is plotted.
#' @param type             Vector indicating type of plot for the statistica power curve
#'                         and the type 1 risk curve.
#' @param pch              Vector indicating plotting symbol for the statistical power curve 
#'                         and the type 1 risk curve.
#' @param lty              Vector indicating line type for the statistical power curve 
#'                         and the type 1 risk curve.
#' @param lwd              Vector indicating line width for the statistical power curve 
#'                         and the type 1 risk curve.                         
#' @param legend           Location of the legend. If \code{FALSE}, legend is omitted.
#' @param bty              Type of box to be drawn around the legend.
#' @param ...              Additional arguments affecting the summary produced.
#' 
#' @author 
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#' Jan Steinfeld \email{jan.steinfeld@@univie.ac.at}
#' 
#' @references
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2009). On designing data-sampling for Rasch model 
#' calibrating an achievement test. \emph{Psychology Science Quarterly, 51}, 370-384.
#'
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2011). A new approach for testing the Rasch model.
#' \emph{Educational Research and Evaluation, 17}, 321-333.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # item parameters
#' ipar2 <- ipar1 <- seq(-3, 3, length.out = 20)
#' # model differential item function (DIF)
#' ipar2[10] <- ipar1[11]
#' ipar2[11] <- ipar1[10]
#' # simulation for b = 100, 200, 300, 400, 500 
#' simres <- pwr.rasch(seq(100, 500, by = 100), ipar = list(ipar1, ipar2))
#' plot(simres)
#' }
plot.pwrrasch <- function(x, plot.sig.level = TRUE, type = c("b", "b"), pch = c(19, 17), 
                          lty = c(1, 3), lwd = c(1, 1), legend = "topleft", bty = "o",  ...) {

  #--------------------------------------------------------------------------------------------------------#
  # Input Check
  
  if (length(x[[1]]) == 1) {
    
    stop("Object pwrrasch contains result of only one b.")
    
  }

  #--------------------------------------------------------------------------------------------------------#
    
  b <- unlist(lapply(x, function(x) x$b))
  
  pwr <- unlist(lapply(x, function(x) x$power))
  
  type1 <- unlist(lapply(x, function(x) x$type1))
  
  sig.level <- unique(unlist(lapply(x, function(x) x$sig.level)))

  ###
  
  plot(b, pwr, type = type[1], lty = lty[1], pch = pch[1], lwd = lwd[1],
       xlab = "", ylab = "", ylim = c(0, 1), axes = FALSE)
  
  axis(1, at = b)    
  axis(2, at = seq(0, 1, by = 0.1))    
  
  mtext("b (number of persons in each group)", side = 1, line = 2.25)
  mtext("Estimated statistical power", side = 2, line = 2.25)
  
  box()
  
  ###
  
  for (i in 1:length(b)) {
    
      lines(c(b[i], b[i]), c(0, pwr[i]), lty = 2, col = "gray70")
      lines(c(b[i], 0), c(pwr[i], pwr[i]), lty = 2, col = "gray70")
    
  }
  
  points(b, pwr, pch = pch[1], type = type[1]) 
  
  
  ###
    
  if (!is.null(type1) & plot.sig.level == TRUE) {    
    
     lines(c(0, max(b)), c(sig.level, sig.level), col = "red2")
    
     points(b, type1, type = type[2], lty = lty[2], pch = pch[2], lwd = lwd[2]) 
    
     if (legend[1] != FALSE) {
       
        if (length(legend) == 1) { 
       
           legend(legend, c("statistical power", "type 1 error risk"), 
                  pch = pch, lty = lty, lwd = lwd, cex = 0.9, bty = bty)
        
        } else {
       
           legend(legend[1], legend[2], c("statistical power", "type 1 error risk"), 
                  pch = pch, lty = lty, lwd = lwd, cex = 0.9, bty = bty)
       
        }
       
     }
     
  } else {
    
    if (legend[1] != FALSE) {
    
       if (length(legend) == 1) { 
        
          legend(legend, "statistical power", 
                 pch = pch[1], lty = lty[1], lwd = lwd[1], cex = 0.9, bty = bty) 
    
      } else {
          
        legend(legend[1], legend[2], "statistical power", 
               pch = pch[1], lty = lty[1], lwd = lwd[1], cex = 0.9, bty = bty)      
      
      }
      
    }
    
  }  
  
}