#' River Chart Scale
#' 
#' This plots the scale of river charts.
#' 
#' 
#' @param length the length of the scale. The length is defined in the same
#' units as the river length. The function will convert this real length into a
#' segment with the same scale as the rivers, and plot it on the river chart.
#' @param label a string defining a scale label indicating the real length the
#' scale segment represents.
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param loc location of scale. One or two values in the range [0, 1] to
#' define left and bottom margin sizes. If \code{loc = NA}, use mouse to locate
#' the arrow. ESC to confirm.
#' @param scl.col scale colour.
#' @param scl.lwd scale line width.
#' @param lbl.cex scale label size.
#' @param lbl.pos scale label position. \code{1} for below, \code{2} for left,
#' \code{3} for above, and \code{4} for right. See \code{par} for details.
#' @param lbl.ofs scale label position offset.
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverDraw}},
#' \code{\link{RiverMap}}, \code{\link{par}}, \code{\link{locator}}.
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
#' RiverScale(2, "2 km", riverlayout, loc = c(0.8, 0.10),lbl.cex = 0.8)
#' 
#' # Use mouse to allocate the river scale
#' ## RiverScale(2, "2 km", riverlayout, lbl.cex = 0.8)
#' 
#' @export RiverScale
RiverScale <- function(length, label, riverlayout,
                       loc = NA, # or loc = c(mar.l, mar.b)
                       scl.col = "black",
                       scl.lwd = 1,
                       lbl.cex = 0.5,
                       lbl.pos = 4,
                       lbl.ofs = 0.5){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # locator
  
  if (all(is.na(loc))){
    loc <- locator()
  }
  
  segments(loc[[1]],loc[[2]], loc[[1]]+length*W.SIZE, loc[[2]], 
           col = scl.col, lty = 1, lwd = scl.lwd)
  segments(loc[[1]], loc[[2]], loc[[1]], loc[[2]] + length*W.SIZE*0.1, 
           col = scl.col, lty = 1, lwd = scl.lwd)
  segments(loc[[1]]+length*W.SIZE, loc[[2]], loc[[1]]+length*W.SIZE, loc[[2]] + length*W.SIZE*0.1, 
           col = scl.col, lty = 1, lwd = scl.lwd)
  
  text(max(loc[[1]],loc[[1]] + length*W.SIZE), loc[[2]], label,
       adj = 0, cex = lbl.cex, pos = lbl.pos, offset = lbl.ofs)
  
}
