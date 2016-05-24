######################################################################
# add.phenos.R
#
# Brian S Yandell
# Ported from http://github.com/byandell/qtlview on 27 apr 2012
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
# Contains: add.rug, qm.approx, myapprox
######################################################################
## This is the only plot routine that refers to same and map.
add.rug <- function(chr, main, maps,
                    p = qm.approx(maps, off.base, chr),
                    use.cM,
                    outer = FALSE,
                    xlim = range(map),
                    bottom.axis = FALSE,
                    side = 1)
{
  bases <- c("cM","Mb")
  base <- bases[2 - use.cM]
  off.base <- bases[1 + use.cM]
  
  ## Add rugs, etc.
  if (length(chr) == 1) {
    ## Get map for chr using proper base.
    map <- maps[[paste(base, "map", sep = ".")]][[chr]]
    ticksize <- ifelse(outer, -0.02, 0.02)

    ## Get plot limits in plotting units.
    usr <- par("usr")

    ## Add grey ticks for non-segregating markers (if available).
    non.seg <- maps[[paste(base, "same", sep = ".")]]
    if(!is.null(non.seg)) {
      non.seg <- non.seg[[chr]]
      rug(non.seg, 0.75 * ticksize, quiet = TRUE, side = side, col = "gray")
      rug(non.seg, 0.75 * ticksize, quiet = TRUE, side = side + 2, col = "gray")
      if(side == 1)
        abline(h = usr[3:4])
      else
        abline(v = usr[1:2])
    }

    rug(map, ticksize, quiet = TRUE, side = side)
    if(bottom.axis) {
      axis(side, pretty(xlim, n = 30), line = ifelse(outer, 0.6, 0))
    }

    rug(map, ticksize, quiet = TRUE, side = side + 2)
    ## This is the culprit.
    axis(side + 2, p$y, p$x, line = ifelse(outer, 0.6, 0))

    usr <- usr[2 * side - c(1,0)]
    usr <- usr[1] - 0.01 * diff(usr[1:2])
    if(use.cM) {
      mtext("cM", side,     1.6, at = usr, adj = 1)
      mtext("Mb", side + 2, 1.6, at = usr, adj = 1)
    }
    else {
      
      mtext("cM", side + 2, 1.6, at = usr, adj = 1)
      mtext("Mb", side,     1.6, at = usr, adj = 1)
    }
    mtext(paste("Chromosome", chr), side, 1.35 + outer)
  }
  title(main, line = 0.5 + 2 * (length(chr) == 1))
}
################################################################
## My approximation routine. Use qm.approx, hide myapprox.
## This cuts off at end of map. Get Karl's approach that extrapolates and add.

qm.approx <- function(maps, base = bases, chr,
                      pos = posn, n.pos = 30,
                      use.qtl = FALSE,
                      ..., non.seg = FALSE)
{
  bases <- c("cM","Mb")
  base <- pmatch(base, bases)[1]
  if(is.na(base))
    stop("base must be cM or Mb")
  
  x <- bases[base]
  y <- bases[-base]
  non.seg <- ifelse(non.seg, "same", "map")
  x <- paste(x, non.seg, sep = ".")
  y <- paste(y, non.seg, sep = ".")

  map.x <- maps[[x]][chr]
  map.y <- maps[[y]][chr]
  posn <- pretty(c(map.x[[1]]), n.pos)
  
  if(use.qtl) {
    ## Need to flesh this out using Aimee's interpolating positions email from 15 nov.
    stop("use.qtl = TRUE is not working yet")
    map.x <- data.frame(...)
    ## interpmap(map.x, map.y)
  }
  else
    myapprox(map.x[[1]], map.y[[1]], pos, ...)
}
################################################################
myapprox <- function(Mb, cM,
  pos = posn, n.pos = 30, ...)
{
  ## Translate Mb to cM within range.
  
  ## Some wierd bug because Mb is of class "A" or "X", but not "numeric".
  posn <- pretty(c(Mb), n.pos)

  ## Adjust pos to be within Mb range.
  tmp <- c(pos)
  tmp <- pmin(max(c(Mb)), pmax(min(c(Mb)), tmp))

  ## Linear interpolation between SNPs.
  p <- approx(c(Mb), c(cM), tmp)

  ## Reset x to be pos.
  p$x <- pos
  p
}
