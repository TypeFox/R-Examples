splooth <- function(x, at, method=c("diff","ratio","wavg")){
  if(is.zoo(x))
    return(xts(splooth(coredata(x), at, method), order.by=index(x)))
  method <- method[1]
  at <- cbind(at, x)[, 1]
  at <- as.logical(at)
  if(!is.matrix(x))
    stop("x must be a matrix.")
  if(ncol(x) < 2)
    stop("x must have at least two columns")
  if(any(n <- is.na(x[which(at), ])))
    stop(paste("x[at, ] cannot have any NA's\n see:", paste(which(n), collapse=",")))
  for(i in 1:(ncol(x) - 1)){
    if(method == "diff"){
      adj <- rep(0, nrow(x))
      adj[at] <- x[at, i] - x[at, i + 1]
      y <- x[, i] - rev(cumsum(rev(adj)))
    }else if(method == "ratio"){
      if(any(x <= 0))
        stop("All x values must be positive if method='ratio'")
      adj <- rep(1, nrow(x))
      adj[at] <- x[at, i] / x[at, i + 1]
      y <- x[, i] * rev(cumprod(rev(adj)))
    }else if(method == "wavg"){
      w <- c(which(at), length(at) + 1)
      tim <- 1:length(at)
      nxt <- w[c(1, (cumsum(at) + 1)[-length(at)])]
      lst <- c(0, cumsum(at)[-length(at)])
      lst[lst != 0] <- w[lst[lst != 0]]
      y <- x[, i] * ((nxt - tim) / (nxt - lst)) + x[, i + 1] * ((tim - lst) / (nxt - lst))
    }else{
      stop("method must be 'diff', 'ratio', or 'wavg'")
    }
    if(i == 1)
      z <- y
    else
      z <- cbind(z, y)
  }
  z <- cbind(z, x[, ncol(x)])
  colnames(z) <- colnames(x)
  rownames(z) <- rownames(x)
  z
}

# this function originally copyright and written by Robert Sams robert@sanctumfi.com 
###############################################################################
# R (http://r-project.org/) Instrument Class Model
#
# Copyright (c) 2009-2012
# Peter Carl, Dirk Eddelbuettel, Jeffrey Ryan, Joshua Ulrich,
# Brian G. Peterson and Garrett See
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: splooth.R 899 2012-01-01 19:00:09Z gsee $
#
###############################################################################
