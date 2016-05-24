#' River Bar-Chart
#' 
#' This plots bars for quantitative data on the river chart.
#' 
#' 
#' @param site a character vector of site names.
#' @param river a vector of rivers on which the sites are located.
#' @param distance a vector. The along-the-river distances between the sites
#' and the mouth of the river.
#' @param value a data frame containing the variables to be shown on the
#' bar-chart.
#' @param riverlayout the output list of \code{RiverLayout} or \code{RiverMap}.
#' @param range bar-chart value range. A vector of two values indicating lower
#' limit and upper limit.
#' @param bar.w relative width of each bar plotted in the diagram. The default
#' value is \code{1}.
#' @param bar.col bar colour.
#' @param bd.col bar border colour.
#' @param pt.shw show location point (\code{TRUE}) or not (\code{FALSE}).
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
#' RiverBar(B.siteaspt$Site, B.siteaspt$River, B.siteaspt$Distance, 
#'          B.siteaspt[4], riverlayout, range = c(0,8), 
#'          bar.col = c("#5381FFFF"), lbl.adj = c(0.5,1.3))
#' 
#' RiverDraw(riverlayout)
#' 
#' RiverBar(B.siteaspt$Site, B.siteaspt$River, B.siteaspt$Distance, 
#'          B.siteaspt[4:5], riverlayout, range = c(0,8), 
#'          bar.col = c("#5381FFFF", "#FF3931FF"), lbl.adj = c(0.5,1.3))
#' 
#' @export RiverBar
RiverBar <- function(site, river, distance, value, riverlayout,
                     range = NA,
                     bar.w = 1,
                     bar.col = NA,
                     bd.col = "black",
                     lbl.cex = 0.7,
                     lbl.adj = c(0.5,2),
                     lbl.ofs = 0.5,
                     lbl.col = "black",
                     lbl.srt = 0,
                     lbl.pos = NULL,
                     lbl.shw = TRUE,
                     pt.shw = FALSE){ 
  
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
  
  # Direction converting
  if (DIRECTION == -1){
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    distance <- length - distance
    X.SITE <- X2[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of sampling sites
  }else{  
    # Site location
    X.SITE <- X1[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of sampling sites
  }
  
  Y.SITE <- Y[match(river, RIVER.DATA$river)] # Y of sampling sites
  
  # Site plot
  if (pt.shw){
    points(X.SITE, Y.SITE) # plot sampling sites
  }
  
  # Bar plotting
  if (all(is.na(range))){
    VALUE.MAX <- max(value)
    VALUE.MIN <- 0
  } else{
    VALUE.MAX <- max(range)
    VALUE.MIN <- min(range)
  }
  
  VALUE.SIZE <- H.SIZE * 0.9/(VALUE.MAX-VALUE.MIN)
  
  # Bar-charts
  W.BAR <- 0.01 * bar.w # bar width, default bar size is 0.01
  value <- data.frame(value)
  N.VALUE <- length(value) # Number of quantitative variables
  
  for (i in 1:N.VALUE){
    rect((X.SITE - W.BAR * N.VALUE/2 + (i-1)*W.BAR), 
         Y.SITE, 
         (X.SITE - W.BAR * N.VALUE/2 + i * W.BAR), 
         (Y.SITE + (value[,i] - VALUE.MIN) * VALUE.SIZE), 
         col = bar.col[i],
         border = bd.col)
  }
  
  if (lbl.shw){
    text(X.SITE, Y.SITE, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
  }
} 