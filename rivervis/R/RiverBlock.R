#' River Block-Chart
#' 
#' This plots blocks to display qualitative data on the river chart.
#' 
#' 
#' @param site a vector of site names.
#' @param river a vector of rivers on which the sites are located.
#' @param distance a vector. The along-the-river distances between the sites
#' and the mouth of the river.
#' @param value a data frame containing the qualitative variables to be shown
#' on the block-chart.
#' @param riverlayout the output list of RiverLayout or RiverMap.
#' @param arrangement a vector indicating the block number for each line.
#' @param pt.shw show location point (\code{TRUE}) or not (\code{FALSE}).
#' @param hw.rat the ratio of block height and width in the plotted diagram.
#' @param h.gap vertical gap size between blocks. By default, the vertical gap
#' is river height * 0.05 in each river chart.
#' @param w.gap horizontal gap size between blocks when there is more than one
#' block in each line. By default, the horizontal gap is largest block width *
#' 0.025.
#' @param block.col block colour.
#' @param block.lwd block line width.
#' @param bd.col block border col.
#' @param par.shw show parameter names (\code{TRUE}) or not (\code{FALSE}).
#' @param par.pos parameter label position. \code{1} for below, \code{2} for
#' left, \code{3} for above, and \code{4} for right. See \code{par} for
#' details.
#' @param par.ofs parameter label position offset.
#' @param par.cex parameter label size.
#' @param par.adj parameter label adjustment. One or two values in [0,1] for x
#' and y (optional) adjustment.
#' @param par.txt parameter name.
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
#' RiverDraw(riverlayout)
#' 
#' RiverBlock(B.sitehm$Site, B.sitehm$River, B.sitehm$Distance, 
#'            B.sitehm[4:9], riverlayout, c(1,1,2,2), 
#'            block.col = fivecolours, lbl.adj = c(0.5,1.3))
#' 
#' RiverDraw(riverlayout)
#' 
#' RiverBlock(B.sitehm$Site, B.sitehm$River, B.sitehm$Distance, 
#'            B.sitehm[4:9], riverlayout, c(1,1,2,2), 
#'            block.col = fivecolours, lbl.adj = c(0.5,1.3),
#'            par.txt = c("ChanVeg", "ChanFlow", "BankVegLeft", 
#'                        "Right", "RipLULeft", "Right"))
#' 
#' 
#' @export RiverBlock
RiverBlock <- function(site, river, distance, value, riverlayout, arrangement,
                       pt.shw = FALSE, # show=1, hide=0
                       hw.rat = 1.5,
                       h.gap = 0.05, # by default, H.GAP is H.SIZE * 0.05
                       w.gap = 0.025, # by default, W.GAP is W.SIZE * 0.025
                       block.col = NA,
                       block.lwd = 1,
                       bd.col = "grey20",                    
                       par.shw = TRUE,
                       par.pos = 2,
                       par.ofs = 1,
                       par.cex =0.6,
                       par.adj = c(1,0.5),
                       par.txt = NA,
                       lbl.shw = TRUE,
                       lbl.cex = 0.7,
                       lbl.adj = c(0.5,2),
                       lbl.ofs = 0.5,
                       lbl.col = "black",
                       lbl.srt = 0,
                       lbl.pos = NULL){
  
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
  
  # Plot block-charts
  N.SITE <- nrow(value) # site/observation number
  N.LINE <- length(arrangement)
  
  H.GAP <- H.SIZE * h.gap # "GAP" is the distance between blocks
  W.GAP <- H.SIZE * w.gap
  
  H.BLOCK <- (H.SIZE-H.GAP)/N.LINE - H.GAP
  # "b" is the height of the rectangle
  W.BLOCK <- H.BLOCK/hw.rat # "a" is the width of the rectangle
  
  Y.BLOCK <- rep(N.LINE:1, arrangement) # repeat the sequence VStr:1, repeating times is provided by VStr
  Y.BLOCK <- H.GAP + (Y.BLOCK-1)*(H.GAP+H.BLOCK)
  X.BLOCK <- sequence(arrangement)-1 # VStr provides the "to"s of the sequence
  
  N.PERLINE <- rep(arrangement, arrangement) # number of blocks per line
  
  for (i in 1:N.SITE){   # draw small blocks site by site
    
    rect((X.SITE[i]-W.BLOCK/2+X.BLOCK*(W.BLOCK+W.GAP)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK),
         (X.SITE[i]-W.BLOCK/2+X.BLOCK*(W.BLOCK+W.GAP)/N.PERLINE+(W.BLOCK-(N.PERLINE-1)*W.GAP)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK + H.BLOCK), 
         col = block.col[as.numeric(value[i,])], border = bd.col, lwd = block.lwd)
  }
  
  # Parameter names
  if (all(is.na(par.txt))){
    PAR.NAME.LIST <- split(colnames(value), rep(1:length(arrangement), arrangement))    
  } else{
    PAR.NAME.LIST <- split(par.txt, rep(1:length(arrangement), arrangement))
  }
  
  PAR.NAMES <- TitlePaste(PAR.NAME.LIST)
  
  X.PAR <- rep(min(X1[RIVER.DATA$row==0],X2[RIVER.DATA$row==0]), N.LINE)
  
  Y.PAR <- Y[RIVER.DATA$row==0] + sort(unique(Y.BLOCK), decreasing = TRUE) + H.BLOCK/2
  
  if (par.shw){
    text(X.PAR, Y.PAR, PAR.NAMES, pos = par.pos, offset = par.ofs, cex = par.cex, adj = par.adj)
  }
  
  if (lbl.shw){
    text(X.SITE, Y.SITE, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
  }
  
  
}
