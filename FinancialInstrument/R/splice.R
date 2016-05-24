splice <- function(x, at){
  if(is.zoo(x))
    return(xts(splice(coredata(x), at), order.by=index(x)))
  x <- as.matrix(x)
  if((nrow(x) %% length(at)) != 0)
    stop("length(at) must be a multiple of nrow(x)")
  at <- as.logical(try(cbind(at, 1:nrow(x))[, 1]))
  if(!any(at))
    return(x)
  w <- which(at)
  n <- max(length(w) - ncol(x) + 1, 1) 
  y <- matrix(NA, nrow=nrow(x), ncol=n)
  for(j in 1:n){
    if(n == j)
      rows <- seq(w[j], nrow(x))
    else
      rows <- seq(w[j], w[j + ncol(x)] - 1)
    cols <- cumsum(at[rows])
    y[rows, j] <- sapply(1:length(rows), function(i){x[rows[i], cols[i]]})
  }
  y
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
# $Id: splice.R 899 2012-01-01 19:00:09Z gsee $
#
###############################################################################
