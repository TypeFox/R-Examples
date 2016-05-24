#' Add an x axis to the current plot
#' 
#' This function adds an x-axis to the current plot. Intended for use when the 
#' plot does not fall along the bottom row, but you plan to put no plot
#' beneath it.
#' 
#' @aliases addxaxis
#@param
#@note %% ~~further notes~~
#' @author Carlisle Rainey (\href{mailto:carlislerainey@@gmail.com}{e-mail},
#' \href{http://www.carlislerainey.com}{website})
#@seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#@references
#' @examples
#' 
#' par(mfrow = c(2,2), mar = c(.75,.75,.75,.75), oma = c(3,3,1,1))
#' eplot(xlim = c(-1, 1), ylim = c(-1, 1), xlab = "X Label")
#' aplot()
#' addxaxis()
#' aplot()
#' 
#' @export addxaxis
#' 
addxaxis <- function() {
  # calculate adjustment factor for axis labels if the plot is a matrix
  deflate <- 1
  if (par("mfg")[3] == 2 & 
        par("mfg")[4] == 2) {
    deflate <- 0.83
  }
  if (par("mfg")[3] > 2 | 
        par("mfg")[4] > 2) {
    deflate <- 0.66
  }
  # add the axis
  axis(side = 1, at = .compactrEnv$plotPar$xat, labels = NA, tck = -.compactrEnv$plotPar$tick.length, 
       lwd = 0, lwd.ticks = 1)
  axis(side = 1, at = .compactrEnv$plotPar$xat, tick = FALSE, line = .compactrEnv$plotPar$xpos, 
       cex.axis =  .9*.compactrEnv$plotPar$text.size)
  mtext(side = 1, text = .compactrEnv$plotPar$xlab, line = .compactrEnv$plotPar$xlabpos, 
        cex = 1*.compactrEnv$plotPar$text.size*deflate)
}