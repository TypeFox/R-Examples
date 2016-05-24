#' Points on River Charts
#' 
#' This function plots scatter points or broken lines on the river chart.
#' 
#' 
#' @param site a vector of site IDs.
#' @param river a vector of river names.
#' @param distance a vector of distances from sites to the river mouth.
#' @param value a vector of values.
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param range point value range. A vector of two values indicating lower
#' limit and upper limit.
#' @param type type of plot. See \code{plot} for details. The default value is
#' "l", which means "lines".
#' @param ln.lwd line width.
#' @param pt.col point or point border colour.
#' @param pt.bg point point background(fill) colour when \code{pt.pch=21:25}.
#' @param pt.pch point style.
#' @param pt.cex point size.
#' @param lbl.cex label size.
#' @param lbl.adj label adjustment. One or two values in [0,1] for x and y
#' (optional) adjustment.
#' @param lbl.ofs label position offset.
#' @param lbl.col label colour.
#' @param lbl.srt label angle.
#' @param lbl.pos label position. \code{1} for below, \code{2} for left,
#' \code{3} for above, and \code{4} for right. See \code{par} for details.
#' @param lbl.shw show labels (\code{TRUE}) or not (\code{FALSE}).
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
#' 
#' RiverDraw(riverlayout)
#' 
#' RiverPoint(B.elevation$Site, B.elevation$River, B.elevation$Distance, 
#'            B.elevation$Elevation, riverlayout)
#' 
#' RiverPoint(B.sitenh4n$Site, B.sitenh4n$River, B.sitenh4n$Distance, 
#'            B.sitenh4n$NH4N_Spring, riverlayout, type = "o", 
#'            pt.col = "#5381FFFF", pt.pch = 21, pt.bg = "lightblue")
#' 
#' RiverPoint(B.sitenh4n$Site, B.sitenh4n$River, B.sitenh4n$Distance, 
#'            B.sitenh4n$NH4N_Autumn, riverlayout, type = "o", 
#'            pt.col = "#FF3931FF", pt.pch = 21, pt.bg = "pink", 
#'            lbl.shw = TRUE)
#' 
#' 
#' @export RiverPoint
RiverPoint <- function(site, river, distance, value, riverlayout,
                       range = NA,
                       type = "l",
                       pt.col = "grey40",
                       pt.bg = "black",
                       pt.pch = 20,
                       pt.cex = 1,
                       lbl.cex = 0.7,
                       lbl.adj = c(0.5,2),
                       lbl.ofs = 0.5,
                       lbl.col = "black",
                       lbl.srt = 0,
                       lbl.pos = NULL,
                       lbl.shw = FALSE,
                       ln.lwd = 1){ 
  
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
  
  # Point plotting 
  if (all(is.na(range))){
    VALUE.MAX <- max(value)
    VALUE.MIN <- 0
  } else{
    VALUE.MAX <- max(range)
    VALUE.MIN <- min(range)
  }
  VALUE.SIZE <- H.SIZE * 0.9/(VALUE.MAX - VALUE.MIN) # a ratio to turn real elevation to plotting scale, assuming the largest value is 0.9*HSIZE
  
  # Direction converting
  if (DIRECTION == -1){
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    distance <- length - distance
    X.VALUE <- X2[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of Elev points on plot
  }else{  
    # Calculate X and Y
    X.VALUE <- X1[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of Elev points on plot
  }
  
  Y.VALUE <- Y[match(river, RIVER.DATA$river)] + (value - VALUE.MIN) * VALUE.SIZE + H.SIZE * 0.1 # use dictionary technique to calculate the Y of Elev points on plot
  
  
  V <- data.frame(river=factor(river), X.VALUE, Y.VALUE)
  
  for (i in RIVER.DATA$river){
    points(V[which(river==i),]$X.VALUE, 
           V[which(river==i),]$Y.VALUE, type=type, col = pt.col, bg = pt.bg, pch = pt.pch, lwd = ln.lwd, cex = pt.cex)
  }
  
  if (lbl.shw){
    X.SITE <- X.VALUE
    Y.SITE <- Y[match(river, RIVER.DATA$river)] # Y of sampling sites
    text(X.SITE, Y.SITE, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
  }
  
}