#' Site of Interest
#' 
#' This plots sites of interest on the river chart.
#' 
#' 
#' @param site a character vector of site names.
#' @param river a vector of rivers on which the sites are located.
#' @param distance a vector. The along-the-river distance between the site and
#' the mouth of the river.
#' @param group a vector. Group names of river locations.
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param pt.pch point style.
#' @param pt.col point border colour.
#' @param pt.bg point background(fill) colour when \code{pt.pch=21:25}.
#' @param pt.cex point size.
#' @param lbl.shw show labels (\code{TRUE}) or not (\code{FALSE}).
#' @param lbl.cex label size.
#' @param lbl.srt label angle.
#' @param lbl.adj label adjustment. One or two values in [0,1] for x and y
#' (optional) adjustment.
#' @param lbl.col label colour.
#' @param lbl.pos label position. \code{1} for below, \code{2} for left,
#' \code{3} for above, and \code{4} for right. See \code{par} for details.
#' @param lbl.ofs label position offset.
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
#' RiverSite(B.soi$SOI, B.soi$River, B.soi$Distance, B.soi$Group, riverlayout, 
#'           pt.bg = c("red","green","yellow"), lbl.shw = FALSE)
#' 
#' RiverDraw(riverlayout)
#' 
#' RiverSite(B.town$Town, B.town$River, B.town$Distance, B.town$Group, 
#'           riverlayout, pt.pch = 22, lbl.shw = FALSE, 
#'           pt.bg = "orange", pt.col = "black")
#' 
#' RiverSite(B.soi$SOI, B.soi$River, B.soi$Distance, B.soi$Group, 
#'           riverlayout, pt.pch = c(25, 24, NA), 
#'           lbl.shw = FALSE, pt.bg = NA, pt.col = "black")
#' 
#' @export RiverSite
RiverSite <- function(site, river, distance, group, riverlayout,
                      pt.pch = 21,
                      pt.col = NA,
                      pt.bg = "red",
                      pt.cex = 1,
                      lbl.cex = 0.5,
                      lbl.srt = 0,
                      lbl.adj = c(0.5,2),
                      lbl.col = "black",
                      lbl.pos = 1,
                      lbl.ofs = 0.5,
                      lbl.shw = TRUE){
  
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
  
  if (DIRECTION == -1){
    X.RIVER <- X2
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    distance <- length - distance
  }else{
    X.RIVER <- X1
  }
  
  
  # Location coordinates
  X.LOC <- X.RIVER[match(river, RIVER.DATA$river)] + distance * W.SIZE # X of locations
  
  Y.LOC <- Y[match(river, RIVER.DATA$river)] # Y of locations
  
  if(length(group)==1){
    points(X.LOC, Y.LOC, type = "p", pch = pt.pch[1], col = pt.col[1], bg = pt.bg[1], cex = pt.cex[1])
    if(lbl.shw){
      text(X.LOC, Y.LOC, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
    }
  }else{
    points(X.LOC, Y.LOC, 
           type = "p", 
           pch = rep(pt.pch, length(group))[c(group)], 
           col = rep(pt.col,length(group))[c(group)], 
           bg = rep(pt.bg,length(group))[c(group)], 
           cex = rep(pt.cex, length(group))[c(group)])
    
    if (lbl.shw){
      text(X.LOC, Y.LOC, labels = site, 
           cex = rep(lbl.cex, length(group))[c(group)], 
           adj = rep(lbl.adj, length(group))[c(group)], 
           col = rep(lbl.col, length(group))[c(group)], 
           srt = lbl.srt,
           pos = lbl.pos,
           offset = lbl.ofs)
    }
  }
  
  
  
  
}
