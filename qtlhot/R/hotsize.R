######################################################################
# hotsize.R
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
# Contains: hotsize, hotsize.scanone, hotsize.highlod,
#           print.hotsize, summary.hotsize, plot.hotsize
######################################################################
hotsize <- function(hotobject, ...) UseMethod("hotsize")

hotsize.scanone <- function(hotobject, lod.thr = NULL, drop.lod = 1.5, ...)
{
  hotsize(highlod(hotobject, lod.thr, drop.lod), lod.thr, ...)
}
hotsize.highlod <- function(hotobject, lod.thr = NULL, window = NULL, quant.level = NULL, ...)
{
  if(length(lod.thr) > 1)
    stop("hotsize only allows one lod.thr value")
  
  scan <- hotobject$chr.pos
  
  highlod <- highlod.thr(hotobject, lod.thr)
  attr(scan, "lod.thr") <- attr(highlod, "lod.thr")
  highlod <- highlod$highlod
  
  if(!nrow(highlod))
    return(NULL)
  
  ## Straight count of LODs above threshold. Includes shoulders and peaks.
  tbl <- table(highlod$row)
  tmp <- rep(0, nrow(scan))
  if(length(tbl))
    tmp[as.numeric(names(tbl))] <- tbl
  scan$max.N <- tmp
  
  ## Smoothed count of peak LODs per chr only above threshold, smoothed by window.
  if(!is.null(window)) {
    window <- round(window)
    if(window > 0) {
      scan$max.N.window <- smooth.neqtl(highlod, scan, window = window)$nqtl
      attr(scan, "window") <- window
    }
  }
  
  if(!is.null(quant.level)) {
    if(is.list(quant.level)) {
      quant.thr <- quant.level$max.N
      lod.thrs <- as.numeric(row.names(quant.thr))
      quant.level <- quant.level$max.lod.quant[,1]
    }
    else
      quant.thr <- NULL
    
    ## If sliding thresholds supplied.
    if(!is.null(lod.thr)) {
      if(!is.null(quant.thr)) {
        m <- which.min(abs(lod.thr - lod.thrs))
        quant.thr <- quant.thr[m,1]
      }
      quant.level <- quant.level[quant.level >= lod.thr & !is.na(quant.level)]
    }
    
    if(length(quant.level)) {
      ## Want to work down quant.level. One way is to see how many are above smallest
      ## and stop if above the index.
      index <- rep(TRUE, nrow(highlod))
      hot.size <- rep(0, nrow(scan))
      for(hot.crit in rev(seq(length(quant.level)))) {
        above <- highlod$lod[index] >= quant.level[hot.crit]
        if(!any(above))
          break;
        tbl <- table(highlod$row[index][above])
        tbl <- tbl[tbl >= hot.crit]
        if(length(tbl)) {
          ## Record hotspot size if above hot.crit for thr, then mask those loci out.
          hot.size[as.numeric(names(tbl))] <- tbl
          index[highlod$row %in% names(tbl)] <- FALSE
        }
        if(!any(index))
          break;
      }
      scan$quant <- hot.size
      attr(scan, "quant.level") <- quant.level
    }
  }
  class(scan) <- c("hotsize", class(scan))
  scan
}
#############################################################################################
print.hotsize <- function(x, ...) print(summary(x, ...))
summary.hotsize <- function(object, ...)
{
  
  cat("hotsize elements: ", paste(names(object)), "\n")

  lod.thr <- attr(object, "lod.thr")
  if(!is.null(lod.thr))
    cat("LOD threshold:", lod.thr, "\n")
  window <- attr(object, "window")
  if(!is.null(window))
    cat("smooth window:", window, "\n")
  quant.level <- attr(object, "quant.level")
  if(!is.null(quant.level)) {
    cat("quantile level summary:\n")
    print(summary(quant.level))
  }
  cat("\n")
  format <- ifelse(ncol(object ==3), "onepheno", "allpeaks")
  keep <- tapply(object$max.N, object$chr, max)
  keep <- object$chr %in% names(keep)[keep > 0]
  if(!any(keep))
    return(invisible())
  object <- object[keep, ]
  NextMethod(object, format = format, ...)
}    
#############################################################################################
plot.hotsize <- function(x, ylab = "counts", quant.axis = pretty(x$max.N),
                         col = c("black","red","blue"), by.chr = FALSE, maps = NULL,
                         title = "", ...)
{
  if(by.chr) {
    for(i in levels(x$chr)) {
      tmp <- x$chr == i
      if(any(tmp))
        plot(x[tmp,], ylab, col = col,
             title = paste("chr", i), ...)

      if(!is.null(maps) & by.chr)
        add.rug(i, "", maps, use.cM = TRUE)
  
    }
  }

  ## Use NextMethod, but repeatedly.
  class(x) <- class(x)[-1]
  ## max.N
  plot(x, lodcolumn = 1, ylab = ylab, col = col[1], ...)
  if(title != "")
    title <- paste("raw (", col[1], ")", sep = "")

  ## max.N.window
  window <- attr(x, "window")
  if(!is.null(window)) {
    plot(x, lodcolumn = 2, col = col[2], add = TRUE, ...)
    if(title != "")
      title <- paste(title, " smoothed(", col[2], ")", sep = "")
  }

  quant.thr <- attr(x, "quant.thr")
  if(!is.null(quant.thr))
    abline(h = quant.thr, lwd = 2, lty = 2, col = col[2])
  
  ## quant
  quant.level <- attr(x, "quant.level")
  if(!is.null(quant.level)) {
    plot(x, lodcolumn = match("quant", names(x)) - 2, col = col[3],
                 add = TRUE, ...)
    if(title != "")
      title <- paste(title, " sliding(", col[3], ")", sep = "")

    ## Add right axis for quantile LOD level.
    if(length(quant.axis)) {
      quant.axis <- pmax(1, quant.axis)
      axis(4, at = quant.axis, labels = round(quant.level[quant.axis], 2), las = 1, cex.axis = 0.9)
      ## mtext("sliding LOD thresholds", 4, 1, cex = 1.5)
    }
  }
  if(title != "")
    mtext(title, 3, 1, cex = 1.5)
  invisible()
}
###################################################################################
smooth.neqtl <- function(highobj, chr.pos, lod.thr = 0, window = 5)
{
  chr <- chr.pos$chr
  pos <- chr.pos$pos
  chr.names <- levels(chr.pos$chr)

  max.hl <- make.maxlod(highobj, chr.pos)
  maxlod.thr.pos <- max.hl$pos
  for(k in chr.names) {
    if(length(max.hl$pos[[k]]))
      maxlod.thr.pos[[k]] <- max.hl$pos[[k]][max.hl$lod[[k]] >= lod.thr]
  }
  
  out <- smoothall(maxlod.thr.pos,thechr = chr.pos$chr, thepos = chr.pos$pos, window = window)
  ## Recover marker information.
  rownames(out) <- rownames(chr.pos)
  out
}
###################################################################################
make.maxlod <- function(highobj, chr.pos)
{
  ## find high LOD and position per chromosome.
  chr.names <- levels(chr.pos$chr)

  tmpfn <- function(x) {
    if(is.null(x))
      0
    else
      max(x, na.rm = TRUE)
  }
  tmpfn2 <- function(x) {
    if(is.null(x))
      NA
    else
      median(x, na.rm=TRUE)
  }
  tmpfn3 <- function(a) {
    is.nan(a) | a==max(a, na.rm=TRUE)
  }

  ## Highlod chromosome. Create it if not supplied.
  hl.chr <- highobj$chr
  if(is.null(hl.chr))
    hl.chr <- chr.pos$chr[highobj$row]
  
  maxlod.hl <- maxlod.pos.hl <- vector("list", length(chr.names))
  names(maxlod.pos.hl) <- chr.names
  for(chr.k in chr.names) {
    ## Subset on chromosome.
    is.chr <- hl.chr == chr.k
    scan.out.bychr <- highobj[is.chr, ]
    
    ## This is kludgey. How to make more efficient?
    scan.out.bychr$phenos <- ordered(scan.out.bychr$phenos, unique(scan.out.bychr$phenos))
    tmp <- unlist(tapply(scan.out.bychr$lod, scan.out.bychr$phenos, tmpfn3))
    if(any(tmp)) {
      scan.out.bychr <- scan.out.bychr[tmp,]
      pos <- chr.pos$pos[highobj$row[is.chr]][tmp]
      
      scan.out.bychr$phenos <- ordered(scan.out.bychr$phenos, unique(scan.out.bychr$phenos))
      ## Find high lod.
      maxlod.hl[[chr.k]] <- tapply(scan.out.bychr$lod, scan.out.bychr$phenos, tmpfn)
      ## Find position of high lod.
      maxlod.pos.hl[[chr.k]] <- tapply(pos, scan.out.bychr$phenos, tmpfn2)
    }
  }
  list(lod = maxlod.hl, pos = maxlod.pos.hl)
}
