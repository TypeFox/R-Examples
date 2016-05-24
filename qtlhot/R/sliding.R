######################################################################
# diag.R
#
# Brian S Yandell
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: slidingbar.create, slidingbar.plot, myplot.err, sliding.bar.plot
######################################################################
slidingbar.create <- function(highobj, quant.level = NULL,
                              restrict.to.levels = FALSE, ...)
{
  if(is.null(quant.level))
    stop("must supply quant.level for sliding bar plot")
  quant.level <- unclass(quant.level)
  
  ## Get matrix of seq(n.quant)
  quant <- quantile(highobj, max.quantile = FALSE)
  n.quant <- nrow(quant)
  rows <- as.numeric(dimnames(quant)[[2]])

  ## Make sure quant.level and quant have same "length".
  l.level <- length(quant.level)
  if(l.level < n.quant) {
    if(restrict.to.levels) {
      n.quant <- l.level
      quant <- quant[seq(n.quant), ]
    }
    else
      quant.level <- c(quant.level,
                       rep(min(quant.level), n.quant - l.level))
  }
  else
    quant.level <- quant.level[seq(n.quant)]
  
  ## set to zero if below quant.level.
  tmpfn <- function(x, q) {
    tmp <- x < q
    if(any(tmp))
      x[tmp] <- 0
    x
  }
  quant[is.na(quant)] <- 0
  quant <- apply(quant, 2, tmpfn, quant.level)

  ## expand quant to whole genome (or at least chr used).
  chrs <- sort(unique(highobj$chr.pos$chr[rows]))
  chr.pos <- highobj$chr.pos[highobj$chr.pos$chr %in% chrs, ]
  expand <- matrix(0, nrow(chr.pos), n.quant)
  expand[rows,] <- t(quant)

  ## probably want to make this into some kind of object?
  out <- data.frame(chr.pos, expand)
  attr(out, "quant.level") <- quant.level

  out
}
slidingbar.plot <- function(x, ...)
{
  x[,-(1:2)] <- 1 * (x[,-(1:2)] > 0)
  col <- c("white","black")
  
  ## Want to borrow from qtlview:::plot.aug.scanone.
  image(seq(nrow(x)), seq(ncol(x) -2), as.matrix(x[,-(1:2)]),
        col = col, xlab = "", ylab = "hotspot size",
        xaxt = "n")

  ## Add chr name
  mtext("Chromosome", 1, 2)
  chr <- levels(x$chr)
  
  n.mar <- table(x$chr)
  wh <- c(0.5, cumsum(n.mar) + 0.5)
  abline(v = wh, xpd = FALSE)
  a <- par("usr")
  abline(v = a[1:2], xpd = FALSE)
  abline(h = a[3:4], xpd = FALSE)
  for (i in 1:length(n.mar))
    axis(side = 1, at = mean(wh[i + c(0, 1)]), labels = chr[i])

  y.axes <- pmax(1, pretty(seq(ncol(x) - 2)))
  quant.level <- round(attr(x, "quant.level"), 2)
  axis(side=4, labels = quant.level[y.axes], at = y.axes, cex.axis = 0.9, las=1)
  invisible()
}
## Generates the sliding bar figure (not an image anymore)
## This is very slow. Rewrite using highlod.
sliding.bar.plot <- function(scan, lod.thr, size.thr, gap=50, y.axes=NULL)
{
  ###
  sliding.bar.matrix <- function(scan, lod.thr, size.thr)
  {
    counts <- count.thr(scan, lod.thr, TRUE)
    out <- matrix(FALSE,nrow(counts),ncol(counts))
    for(i in order(size.thr))
      out[ i, counts[i,] > size.thr[i] ] <- TRUE
    out
  }
  ###
  matrix.to.plot <- function(M, map, myrug)
  {
    pM <- matrix(FALSE, nrow(M), ceiling(max(myrug)))
    rrug <- as.numeric(round(myrug,0))
    for(i in 1:ncol(M))
      pM[,rrug[i]] <- M[,i]
    pM
  }
  ###
  create.rug <- function(map, gap)
  {
    chrs <- unique(map[,1])
    nchrs <- length(chrs)
    maxpos <- rep(NA,nchrs)
    myrug <- map[map[,1]==1,2]
    chr.legend.pos <- median(myrug)
    for(i in 2:length(chrs)){
      aux <- max(myrug) + gap
      myrug <- c(myrug, map[map[,1]==i,2] + aux)
      chr.legend.pos <- c(chr.legend.pos, median(map[map[,1]==i,2] + aux))
    }
    list(myrug=myrug, chr.legend.pos=chr.legend.pos)
  }
  ###
  N <- max(size.thr)
  sbm <- sliding.bar.matrix(scan, lod.thr, size.thr)
  map <- scan[,1:2]
  myrug <- create.rug(map, gap)
  M <- matrix.to.plot(sbm, map, myrug[[1]])
  lod.thr <- as.character(round(lod.thr,2))
  xaxis <- c(1:ncol(M))
  par(mar=c(5, 4, 4, 5) + 0.1) 
  plot(xaxis, xaxis, type="n", ylim=c(0,N), xaxt="n", xlab="Chromosome",
       yaxt="n", ylab="Hotspot size", cex.lab=1.5)
  rug(myrug[[1]], 0.02, quiet=TRUE)
  axis(side=1, labels=as.character(unique(map[,1])), at=myrug[[2]], 
       cex.axis=1.5, tick=FALSE)
  if(is.null(y.axes))
    y.axes <- quantile(c(1:N), c(0, 0.25, 0.50, 0.75, 1))
  axis(side=2, labels=as.character(y.axes), at=y.axes, cex.axis=0.9, las=1)
  axis(side=4, labels=lod.thr[y.axes], at=y.axes, cex.axis=0.9, las=1)
  mtext("LOD threshold",side=4,cex=1.4,line=3.5,adj=0.545)
  ## This is the slow part. Better to do as image?
  for(i in 1:N){
    for(j in 1:ncol(M)){
      if(M[i,j]) segments(x0=j, x1=j, y0=i-1/2, y1=i+1/2, lwd=0.1)
    }
  }
}

