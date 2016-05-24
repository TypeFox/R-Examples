#' River Chart Plotting
#' 
#' This plots the river chart according to the output list provided by
#' \code{RiverLayout}.
#' 
#' 
#' @param riverlayout the output list of \code{RiverLayout}.
#' @param bd.col border colour.
#' @param bg.col background colour.
#' @param ln.col lead line colour.
#' @param ln.lty lead line style.
#' @param ln.lwd lead line width.
#' @param pt.shw show anchor point (\code{TRUE}) or not (\code{FALSE}). Anchor
#' points represent the locations of the river mouths.
#' @param pt.col anchor point colour.
#' @param pt.pch anchor point style.
#' @param pt.bg anchor point background(fill) colour when \code{pch=21:25}.
#' @param pt.cex anchor point size.
#' @param pt.lwd anchor point border width.
#' @param mar.t top margin size.
#' @param mar.b bottom margin size.
#' @param mar.l left margin size.
#' @param mar.r right margin size.
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverMap}},
#' \code{\link{par}}.
#' @keywords hplot
#' @examples
#' 
#' 
#' data(Ballinderry)
#' 
#' # River flows right
#' riverlayout <- RiverLayout(B.river$River,B.river$Length,
#'                            B.river$Parent,B.river$Position,
#'                            B.river$Distance, direction = -1)
#' RiverDraw(riverlayout)
#' 
#' # River flows left
#' riverlayout.left <- RiverLayout(B.river$River,B.river$Length,
#'                                 B.river$Parent,B.river$Position,
#'                                 B.river$Distance)
#' 
#' RiverDraw(riverlayout.left)
#' 
#' @export RiverDraw
RiverDraw <- function(riverlayout,
                      bd.col = "black",
                      ln.col = "grey40",
                      ln.lty = 3,
                      ln.lwd = 1,
                      bg.col = "grey80",
                      pt.shw = TRUE,
                      pt.col = "black",
                      pt.pch = 20,
                      pt.bg = "black",
                      pt.cex = 1,
                      pt.lwd = 1,
                      mar.t = 0.05,  
                      mar.b = 0.05,
                      mar.l = 0.2,
                      mar.r = 0.1
){
  
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
  
  
  par(mar=c(0,0,0,0))
  
  # Plot new sheet
  plot.new()
  
  # plotting margin
  par(usr = c(-mar.l, 1+mar.r, -mar.b, 1+mar.t))
  
  # Plot river rectangles
  rect(X1, Y, X2, Y+H.SIZE, col = bg.col, border = bd.col) # draw river rectangles
  
  # Plot lead line
  Y.PARENT <- Y[match(RIVER.DATA$parent, RIVER.DATA$river)] # Y of Parent of each river, using dictionary technique
  
  segments(X1, Y, X1, Y.PARENT, col = ln.col, lty = ln.lty, lwd = ln.lwd)
  
  # Plot river rectangle-frames
  rect(X1, Y, X2, Y+H.SIZE, col = NA, border = bd.col) # draw the frame of river rectangles, in case they have been covered by leadlines
  
  # Plot anchor points
  if (pt.shw){
    points(X1, Y.PARENT, pch=pt.pch, col = pt.col, bg = pt.bg, cex = pt.cex, lwd = pt.lwd) # plot anchor points
  }
  
}
