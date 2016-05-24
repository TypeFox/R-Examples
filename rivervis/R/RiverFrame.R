#' River Frame Plotting
#' 
#' This plots river frames, lead lines and archor points.
#' 
#' 
#' @param riverlayout the output list of \code{RiverLayout}.
#' @param bd.col border colour.
#' @param bd.shw show boarders (\code{TRUE}) or not (\code{FALSE}).
#' @param ln.col lead line colour.
#' @param ln.shw show lead lines (\code{TRUE}) or not (\code{FALSE}).
#' @param ln.lty lead line style.
#' @param ln.lwd lead line width.
#' @param pt.shw show anchor points (\code{TRUE}) or not (\code{FALSE}). Anchor
#' points represent the locations of the river mouths.
#' @param pt.col anchor point colour.
#' @param pt.pch anchor point style.
#' @param pt.bg anchor point background(fill) colour when \code{pch=21:25}.
#' @param pt.cex anchor point size.
#' @param pt.lwd anchor point border width.
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverMap}},
#' \code{\link{par}}.
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
#' RiverFrame(riverlayout, bd.col = "green", 
#'            pt.col = "red", ln.col = "orange")
#' 
#' @export RiverFrame
RiverFrame <- function(riverlayout,
                       ln.shw = T,
                       ln.col = "grey40",
                       ln.lty = 3,
                       ln.lwd = 1,
                       pt.shw = T,
                       pt.col = "black",
                       pt.pch = 20,
                       pt.bg = "black",
                       pt.cex = 1,
                       pt.lwd = 1,
                       bd.shw = T,
                       bd.col = "black"){
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
  
  # Plot frame
  if (bd.shw){
    rect(X1, Y, X2, Y+H.SIZE, border = bd.col) # draw river rectangles    
  }
  
  # Plot lead line
  Y.PARENT <- Y[match(RIVER.DATA$parent, RIVER.DATA$river)] # Y of Parent of each river, using dictionary technique
  
  if (ln.shw){
    segments(X1, Y, X1, Y.PARENT, col = ln.col, lty = ln.lty, lwd = ln.lwd)    
  }
  
  # Plot anchor points
  if (pt.shw){
    points(X1, Y.PARENT, pch=pt.pch, col = pt.col, bg = pt.bg, cex = pt.cex, lwd = pt.lwd) # plot anchor points
  }
}
