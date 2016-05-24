#' River Block Chart Plotting
#' 
#' This function plots a river block chart to display qualitative data without
#' the topological structure of the river network. The function does not
#' require the output list from \code{RiverLayout} or \code{RiverMap}.
#' 
#' 
#' @param site a character vector of site names.
#' @param river a vector of rivers to which the sites belong.
#' @param distance a vector. The along-the-river distance between the site and
#' the mouth of the river.
#' @param value a data frame containing the qualitative variables to be shown
#' on the river block-chart.
#' @param arrangement a vector indicating the block number for each line.
#' @param h.gap vertical gap size between blocks in the plot. By default, the
#' vertical gap is block height * 0.1
#' @param w.gap horizontal gap size between sites. By default, the horizontal
#' gap is largest block width * 0.25.
#' @param w.gap.s horizontal gap size between small blocks when there are more
#' than one block in each line. By default, the horizontal gap is largest block
#' width * 0.1.
#' @param r.gap horizontal gap size between rivers. By default, gap between
#' rivers is block width * 0.25.
#' @param block.col a vector of block colours. The length of vector should be
#' as the same as the the number of levels of \code{value}.
#' @param block.lwd a value of block line width.
#' @param border.col a value of block border colour.
#' @param bg.col a value of river background colour.
#' @param mar a value of smallest margin size.
#' @param hw.rat the ratio of block height and block width.
#' @param site.shw show site names (\code{TRUE}) or not (\code{FALSE}).
#' @param site.pos site position. \code{1} for below, \code{2} for left,
#' \code{3} for above, and \code{4} for right. See \code{par} for details.
#' @param site.ofs site position offset.
#' @param site.cex site name size.
#' @param site.col site colour.
#' @param site.order order of sites within each river. Alphabetical order
#' (\code{"A"}), river flow left (\code{"L"}), river flow right (\code{"R"}).
#' @param site.srt site label rotation in degrees.
#' @param rvr.shw show river labels (\code{TRUE}) or not (\code{FALSE}).
#' @param rvr.t.b location of river label. \code{"t"} for top (above) and
#' \code{"b"} for below (bottom).
#' @param rvr.ofs river label position offset.
#' @param rvr.cex river label size.
#' @param rvr.col river label colour.
#' @param rvr.order order of rivers. Alphabetical order (\code{NA}) or a vector
#' of custom order.
#' @param rvr.srt river label rotation in degrees.
#' @param par.shw show parameter labels (\code{TRUE}) or not (\code{FALSE}).
#' @param par.pos parameter label position. \code{1} for below, \code{2} for
#' left, \code{3} for above, and \code{4} for right. See \code{par} for
#' details. It overrides \code{par.adj} if par.pos is not \code{NULL}.
#' @param par.ofs parameter label position offset.
#' @param par.cex parameter label size.
#' @param par.col parameter label colour.
#' @param par.adj parameter label adjustment. One or two values in [0,1] for x
#' and y (optional) adjustment.
#' @param par.txt parameter name.
#' @author Feng Mao
#' @seealso \code{\link{par}}.
#' @keywords hplot
#' @examples
#' 
#' 
#' data(Ballinderry)
#' 
#' RiverBlockChart(B.sitehm$Site, B.sitehm$River, B.sitehm$Distance, 
#'                 B.sitehm[4:9],  c(1,1,2,2), mar = 0.15, site.ofs = 1, 
#'                 site.cex = 0.7, site.order = "R",
#'                 block.col = fivecolours)
#' 
#' RiverBlockChart(B.sitehm$Site, B.sitehm$River, B.sitehm$Distance, 
#'                 B.sitehm[4:9],  c(1,1,2,2), mar = 0.15, 
#'                 site.ofs = 1, site.cex = 0.7, 
#'                 rvr.order = c("Rock", "Killymoon-Claggan", "Ballinderry", 
#'                               "Ballymully", "Kildress", "Kingsmill", 
#'                               "Lissan", "Tulnacross"),
#'                 block.col = fivecolours)
#' 
#' RiverBlockChart(B.sitehm$Site, B.sitehm$River, B.sitehm$Distance, 
#'                 B.sitehm[4:9],  c(1,1,2,2), mar = 0.15, site.ofs = 1, 
#'                 site.cex = 0.7, site.order = "R",
#'                 par.txt = c("ChanVeg", "ChanFlow", "BankVegLeft", 
#'                             "Right", "RipLULeft", "Right"),
#'                 block.col = fivecolours)
#' 
#' @export RiverBlockChart
RiverBlockChart <- function(site, river, distance, value, arrangement,
                            h.gap = 0.1, # by default, H.GAP is H.BLOCK * 0.5
                            w.gap = 0.25, # by default, W.GAP is W.BLOCk * 0.25
                            w.gap.s = 0.1, # by default, gap between small blocks is W.BLOCK *0.1
                            r.gap = 0.25, # by default, gap between rivers is W.BLOCK * 0.25
                            block.col = NA,
                            block.lwd = 1,
                            border.col = "grey20",
                            bg.col = "lightgrey",
                            mar = 0.1, # smallest margin
                            hw.rat = 1.5,
                            site.shw = TRUE,
                            site.pos = 1,
                            site.ofs = 1.5,
                            site.cex = 0.5,
                            site.col = "black",
                            site.order = "A", # alphabetical order ("A"), river flow left ("L"), river flow right ("R")
                            site.srt = 0,
                            rvr.shw = TRUE,
                            rvr.ofs = 1.5,
                            rvr.cex = 0.7,
                            rvr.col = "black",
                            rvr.t.b = "b",
                            rvr.order = NA, # alphabetical order (NA) or custom order (a vector)
                            rvr.srt = 0,
                            par.shw = TRUE,
                            par.pos = 2,
                            par.ofs = 1,
                            par.cex =0.6,
                            par.adj = c(1,0.5),
                            par.col = "black",
                            par.txt = NA){ 
  
  # plot new
  par(mar=c(0,0,0,0))
  plot.new()
  
  # layout
  N.SITE <- nrow(value) # site/observation number
  N.LINE <- length(arrangement)
  N.RIVER <- length(unique(river))
  
  if((N.SITE+N.RIVER-1)>N.LINE){
    MARGIN.LEFT <- mar
    MARGIN.RIGHT <- 1- mar
    WIDTH <- 1-mar*2
    W.BLOCK <- WIDTH/((1+w.gap)*N.SITE + N.RIVER*w.gap + (N.RIVER-1)*r.gap)
    H.BLOCK <- W.BLOCK * hw.rat
    HEIGHT <- H.BLOCK * (N.LINE + h.gap * N.LINE +h.gap)
    MARGIN.TOP <- (1-HEIGHT)/2
    MARGIN.BOTTOM <- MARGIN.TOP
  }else{
    MARGIN.TOP <- 1- mar
    MARGIN.BOTTOM <- mar
    HEIGHT <- 1-mar*2
    H.BLOCK <- HEIGHT/(N.LINE + h.gap * N.LINE +h.gap)
    W.BLOCK <- H.BLOCK/hw.rat
    WIDTH <- (W.BLOCK + w.gap*W.BLOCK) * N.SITE + N.RIVER*w.gap*W.BLOCK + (N.RIVER-1)*r.gap*W.BLOCK
    MARGIN.LEFT <- (1-WIDTH)/2
    MARGIN.RIGHT <- 1- MARGIN.LEFT
    
  }
  
  H.GAP <- H.BLOCK * h.gap # "GAP" is the distance between blocks
  W.GAP <- W.BLOCK * w.gap
  W.GAP.S <- W.BLOCK * w.gap.s
  R.GAP <- W.BLOCK * r.gap
  
  
  if(all(is.na(rvr.order))){
    rvr.order <- sort(unique(river))
  }
  
  # Site coordinates
  
  LOC.RIVER.CAL <- data.frame(river = sort(unique(river)), number = matrix(table(river)))
  
  LOC.RIVER.CAL$river <- factor(LOC.RIVER.CAL$river, levels = rvr.order)
  
  LOC.RIVER.CAL <- LOC.RIVER.CAL[order(LOC.RIVER.CAL$river),]
  
  for(i in N.RIVER:1){
    LOC.RIVER.CAL$accnumber[i] <- sum(LOC.RIVER.CAL$number[1:i-1])
  }
  
  LOC.RIVER.CAL$gapnum <- 0:(N.RIVER-1)
  LOC.RIVER.CAL$x <- LOC.RIVER.CAL$gapnum*(W.GAP + R.GAP) + LOC.RIVER.CAL$accnumber * (W.BLOCK + W.GAP)
  
  A <- data.frame(site, river, distance, NumberinGroup = rep(0,N.SITE), value)
  A$river <- factor(A$river, levels = rvr.order)
  
  if (site.order == "A"){
    A <- A[order(A$river, A$site),]
  }else if (site.order == "L"){
    A <- A[order(A$river, A$distance),]
  }else if (site.order == "R"){
    A <- A[order(A$river, -A$distance),]
  }else{
    
  }
  
  for (i in river){
    k = 1
    for (j in 1:N.SITE){
      if (A$river[j]==i){
        A[j, "NumberinGroup"] <- k
        k <- k + 1
      }
    }
  }
  
  X.RIVER <- LOC.RIVER.CAL$x[match(A$river,LOC.RIVER.CAL$river)]
  
  NumberinGroup <- A$NumberinGroup
  
  
  X.SITE <-  MARGIN.LEFT + X.RIVER + (NumberinGroup-1)*(W.BLOCK + W.GAP) + W.GAP
  
  Y.SITE <- rep((H.GAP + MARGIN.BOTTOM), N.SITE) # Y of sampling sites
  
  
  # Background
  
  X1.RIVER.BG <- MARGIN.LEFT + LOC.RIVER.CAL$x
  Y1.RIVER.BG <- MARGIN.BOTTOM + rep(0,N.RIVER) - 3 * H.GAP
  X2.RIVER.BG <- MARGIN.LEFT + LOC.RIVER.CAL$x + LOC.RIVER.CAL$number*(W.BLOCK+W.GAP)+W.GAP
  Y2.RIVER.BG <- MARGIN.BOTTOM +  rep(HEIGHT,N.RIVER) + 3 * H.GAP
  
  rect(X1.RIVER.BG, Y1.RIVER.BG, X2.RIVER.BG, Y2.RIVER.BG, border = NA, col = bg.col)
  
  
  # Block coordinates
  Y.BLOCK <- rep(N.LINE:1, arrangement) # repeat the sequence VStr:1, repeating times is provided by VStr
  Y.BLOCK <- (Y.BLOCK-1)*(H.GAP+H.BLOCK)
  X.BLOCK <- sequence(arrangement)-1 # VStr provides the "to"s of the sequence
  
  N.PERLINE <- rep(arrangement, arrangement) # number of blocks per line
  
  VALUE.BLOCK <- A[(ncol(A)-ncol(value)+1):ncol(A)]
  
  for (i in 1:N.SITE){   # draw small blocks site by site
    
    rect((X.SITE[i]+X.BLOCK*(W.BLOCK+W.GAP.S)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK),
         (X.SITE[i]+X.BLOCK*(W.BLOCK+W.GAP.S)/N.PERLINE+(W.BLOCK-(N.PERLINE-1)*W.GAP.S)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK + H.BLOCK), 
         col = block.col[as.numeric(VALUE.BLOCK[i,])], border = border.col, lwd = block.lwd)
    
  }
  
  # Site name
  if (site.shw){
    text(X.SITE+W.BLOCK/2, Y.SITE, A$site, pos= site.pos, offset = site.ofs, cex = site.cex, col = site.col, srt = site.srt)
  }
  
  # River name
  if (rvr.shw){
    if(!all(is.na(rvr.order))){
      RIVER.NAME <- rvr.order
    } else{
      RIVER.NAME <- sort(unique(river))
    }
    
    if (rvr.t.b == "b"){
      text((X1.RIVER.BG + X2.RIVER.BG)/2, Y1.RIVER.BG, RIVER.NAME, pos = 1, offset = rvr.ofs, cex = rvr.cex, col = rvr.col, srt = rvr.srt)
    } 
    if (rvr.t.b == "t"){
      text((X1.RIVER.BG + X2.RIVER.BG)/2, Y2.RIVER.BG, RIVER.NAME, pos = 3, offset = rvr.ofs, cex = rvr.cex, col = rvr.col, srt = rvr.srt)
    }
    
  }
  
  # Parameter names
  if (all(is.na(par.txt))){
    PAR.NAME.LIST <- split(colnames(value), rep(1:length(arrangement), arrangement))    
  } else{
    PAR.NAME.LIST <- split(par.txt, rep(1:length(arrangement), arrangement))
  }
  
  PAR.NAMES <- TitlePaste(PAR.NAME.LIST)
  
  X.PAR <- rep(min(X.SITE), N.LINE)
  
  Y.PAR <- sort(unique(min(Y.SITE)+Y.BLOCK + H.BLOCK/2), decreasing = TRUE)
  
  if (par.shw){
    text(X.PAR, Y.PAR, PAR.NAMES, pos = par.pos, offset = par.ofs, cex = par.cex, adj = par.adj, col = par.col)
  }
  
  
}
