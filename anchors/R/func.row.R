## File:   func.row.R
## Author: Jonathan Wand
## Date:   2002-
##
## Created :  2002-06-27
## Modified:  $Date: 2005/08/10 19:06:50 $
## Revision:  $Revision: 1.1 $
## RCS-ID:    $Id: func.row.R,v 1.1 2005/08/10 19:06:50 jwand Exp $
##
## Purpose:
##   Cumulative functions on rows
##   Faster than using apply(,,cumsum)
##
row.cumsum <- function(a) {
  nrow <- dim(a)[1]
  ncol <- dim(a)[2]
  b <- matrix(NA,nrow=nrow,ncol=ncol)
  for (i in 1:nrow)
    b[i,] <- cumsum(a[i,])
  return(b)
}
row.cummax <- function(a) {
  nrow <- dim(a)[1]
  ncol <- dim(a)[2]
  b <- matrix(NA,nrow=nrow,ncol=ncol)
  for (i in 1:nrow)
    b[i,] <- cummax(a[i,])
  return(b)
}

