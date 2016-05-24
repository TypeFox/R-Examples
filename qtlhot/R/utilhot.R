######################################################################
# util.R
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
# Contains: lininterp, findpeaks
######################################################################

## see ~/p/private/diabetes1/diabetes10/scan.perm/func1.R
lininterp <- function(x, y, xnew, ynew)
{
  if(!missing(ynew)) {
    if(!missing(xnew)) 
      stop("Give one of xnew or ynew, not both.")
    return(lininterp(y, x, ynew))
  }

  if(missing(xnew))
    stop("Give one of xnew or ynew.")

  wh <- is.nan(x) | is.nan(y)
  if(any(wh)) {
    x <- x[!wh]
    y <- y[!wh]
  }
  if(length(x) != length(y))
    stop("x and y must have the same length")
  if(length(x) < 2)
    stop("x and y must have length > 1")
  if(any(diff(x) < 0) || any(diff(y) < 0))
    stop("x and y must be non-decreasing")
  
  xd <- diff(range(x))
  yd <- diff(range(y))

  z <- rep(NA, length(xnew))
  seen <- rep(FALSE, length(xnew))

  # on top
  wh <- match(xnew, x)
  if(any(!is.nan(wh))) {
    z[!is.nan(wh)] <- y[wh[!is.nan(wh)]]
    seen[!is.nan(wh)] <- TRUE
  }
    
  # to the right
  wh <- !is.nan(xnew) & xnew > max(x) & !seen
  if(any(wh)) {
    z[wh] <- (xnew[wh]-max(x))/xd*yd+max(y)
    seen[wh] <- TRUE
  }
  
  # to the left
  wh <- !is.nan(xnew) & xnew < min(x) & !seen
  if(any(wh)) {
    z[wh] <- min(y) - (min(x)-xnew[wh])/xd*yd
    seen[wh] <- TRUE
  }

  # in between
  wh <- !is.nan(xnew) & xnew >= min(x) & xnew <= max(x) & !seen
  if(any(wh)) {
    seen[wh] <- TRUE

    for(i in which(wh)) {
      xl <- max(x[x <= xnew[i]])
      xr <- min(x[x >= xnew[i]])
      yl <- max(y[x <= xnew[i]])
      yr <- min(y[x >= xnew[i]])
      z[i] <- (xnew[i]-xl)/(xr-xl)*(yr-yl) + yl
    }
  }
  
  z
}

findpeaks <- function(results, lodcolumn=7, window=5, n.peaks=20)
{
  results$chr <- as.character(results$chr)
  for(i in 1:n.peaks) {
    temp <- max(results, lodcolumn=lodcolumn)
    if(i==1) 
      output <- temp
    else
      output <- rbind(output, temp)
    
    chr <- as.character(temp[[1]])
    pos <- temp[[2]]
    results[results$chr==chr & results$pos >= pos-window & results$pos <= pos+window,lodcolumn+2] <- NA
  }
  output
}
      
######################################################################
get.chr.pos <- function(cross)
{
  if(is.null(cross$geno[[1]]$prob))
    stop("must first run calc.genoprob with proper settings")
  ncross <- cross
  ncross$pheno <- data.frame(trait = rnorm(nind(cross)))
  scanone(ncross, pheno.col= find.pheno(ncross, "trait"))[,1:2]
}
