#' River Reach Plotting
#' 
#' This highlights river reaches on the river chart.
#' 
#' 
#' @param reach a vector of reach names.
#' @param river a vector of rivers to which the reaches belong.
#' @param from a numeric vector of starting points.
#' @param to a numeric vector of ending points.
#' @param group a vector of reach group names. This indicates to which group
#' the reaches belong.
#' @param style a vector of reach styles. The value of "style" denotes the
#' location of reach lines. Especially, \code{0} denotes "on axis" and
#' \code{99} means "the reach is presented as a band rather than a line".
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param rea.pos a vector of absolute positions of lines. The values range in
#' [0,1].
#' @param rea.col line colour.
#' @param rea.lty line style.
#' @param rea.lwd line width.
#' @param rea.den the density of shading lines, in lines per inch. See
#' \code{rect}.
#' @param bd.col colour of river chart frames.
#' @param ln.shw show lead lines (\code{TRUE}) or not (\code{FALSE}).
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
#' @author Feng Mao
#' @seealso \code{\link{RiverLayout}}, \code{\link{RiverDraw}},
#' \code{\link{RiverMap}}, \code{\link{par}}.
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
#' RiverReach(B.reach$Reach, B.reach$River, B.reach$From, 
#'            B.reach$To, B.reach$Group, B.reach$Style, riverlayout,
#'            rea.lwd = 5)
#' 
#' RiverReach(B.reach$Reach, B.reach$River, B.reach$From, 
#'            B.reach$To, B.reach$Group, 2, riverlayout,
#'            rea.col = "darkred", rea.lwd = 5)
#' 
#' @export RiverReach
RiverReach <- function(reach, river, from, to, group, style, riverlayout,
                       rea.pos = NA,       # absolute positions of lines
                       rea.col = "lightblue",
                       rea.lty = 1,
                       rea.lwd = 1,
                       rea.den = NULL,
                       bd.col = "black",
                       ln.shw = T,
                       ln.col = "grey40",
                       ln.lty = 3,
                       ln.lwd = 1,
                       pt.shw = T,
                       pt.col = "black",
                       pt.pch = 20,
                       pt.bg = "black",
                       pt.cex = 1,
                       pt.lwd = 1){
  
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
  
  # Lines or Bands
  Style.Line <- style
  Style.Line[Style.Line == 99] <- NA
  
  Style.Band <- style
  Style.Band[Style.Band != 99] <- NA
  
  if (length(style[style<99]) > 0){
    N.LINE <- max(style[style<99])
  }else{
    N.LINE <- 0
  }
  N.OBS <- nrow(style)
  
  
  # Direction converting
  
  if (DIRECTION == -1){
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    from <- length - from
    to <- length - to
    X.FROM <- X2[match(river, RIVER.DATA$river)] + from * W.SIZE # X of reach
    X.TO <- X2[match(river, RIVER.DATA$river)] + to * W.SIZE # X of reach
  }else{  
    # Calculate X and Y
    X.FROM <- X1[match(river, RIVER.DATA$river)] + from * W.SIZE # X of reach
    X.TO <- X1[match(river, RIVER.DATA$river)] + to * W.SIZE # X of reach
    
  }
  
  Y.REACH.BASE <- Y[match(river, RIVER.DATA$river)] # Y base of locations
  
  Y.REACH.BAND <- Y.REACH.BASE + Style.Band-99
  
  # draw bands and lines
  
  rect(X.FROM, Y.REACH.BAND, X.TO, (Y.REACH.BAND + H.SIZE),
       border = bd.col, density = rea.den,
       col = rep(rea.col, nlevels(group))[c(group)]
  )
  
  if (all(is.na(rea.pos))){
    
    Y.REACH.LINE <- Y.REACH.BASE + H.SIZE/(N.LINE+1)*Style.Line
    
    
  }else{
    
    rea.pos[which(style == 99)] <- NA
    
    Y.REACH.LINE <- Y.REACH.BASE + rea.pos * H.SIZE
    
  }
  
  segments(X.FROM, Y.REACH.LINE, X.TO, Y.REACH.LINE, 
           col = rep(rea.col, nlevels(group))[c(group)],
           lty = rep(rea.lty, nlevels(group))[c(group)], 
           lwd = rep(rea.lwd,nlevels(group))[c(group)])
  
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
