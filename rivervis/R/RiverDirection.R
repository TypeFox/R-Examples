#' River Direction Arrow
#' 
#' This plots river flow direction arrow on river charts.
#' 
#' 
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param loc location of arrow. One or two values in the range [0, 1] for left
#' and bottom margin sizes. If \code{loc = NA}, use mouse to locate the arrow.
#' ESC to confirm.
#' @param arw.length arrow length.
#' @param arw.lty arrow line style.
#' @param arw.lwd arrow line width.
#' @param arw.angle arrow head angle.
#' @param arw.col arrow colour.
#' @param label label of the arrow.
#' @param lbl.cex label size.
#' @param lbl.pos label position.
#' @param lbl.ofs label position offset.
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverDraw}},
#' \code{\link{RiverMap}}, \code{\link{par}}, \code{\link{locator}},
#' \code{\link{arrows}}.
#' @keywords hplot
#' @examples
#' 
#' 
#' data(Ballinderry)
#' 
#' riverlayout <- RiverLayout(B.river$River,B.river$Length,
#'                            B.river$Parent,B.river$Position,
#'                            B.river$Distance, direction = -1)
#' RiverDraw(riverlayout)
#' 
#' RiverDirection(riverlayout, arw.length = 0.03, 
#'                loc = c(0.8, 0.05), lbl.cex = 0.8)
#' 
#' # Use mouse to allocate the flow direction sign
#' ## RiverDirection(riverlayout, arw.length = 0.03, lbl.cex = 0.8)
#' 
#' 
#' @export RiverDirection
RiverDirection <- function(riverlayout, 
                           loc = NA, # or loc = c(mar.l, mar.b)
                           arw.length = 0.05,
                           arw.lty = 1,
                           arw.lwd = 1,
                           arw.angle = 30,
                           arw.col = "black",
                           label = "Flow direction",
                           lbl.cex = 0.5,
                           lbl.pos = 4,
                           lbl.ofs = 0.5){
  
  # Data transfer
  
  DIRECTION <- riverlayout[[9]]
  
  # locator
  
  if (all(is.na(loc))){
    loc <- locator()
  }
  
  arrows(loc[[1]], loc[[2]], loc[[1]] + arw.length, loc[[2]], 
         length = arw.length * 2, code = 1.5-DIRECTION/2, lty = arw.lty, lwd = arw.lwd, angle = arw.angle, col = arw.col)
  
  text(max(loc[[1]],loc[[1]] + arw.length), loc[[2]], label,
       adj = 0, cex = lbl.cex, pos = lbl.pos, offset = lbl.ofs)
}
