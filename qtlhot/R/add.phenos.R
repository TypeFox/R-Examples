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
# Contains: add.phenos
######################################################################
add.phenos <- function(cross, newdata = NULL, index = NULL)
{
  if(!is.null(newdata)) {
    if(any(names(newdata) %in% names(cross$pheno)))
      warning("some cross phenotypes overwritten with new data")

    if(is.null(index)) {
      n.ind <- nind(cross)
      if(nrow(newdata) != n.ind)
        stop(paste("newdata must have number of rows (",
                   nrow(newdata), ") as cross individuals (",
                   n.ind, ")", sep = ""))

      ## Must assume newdata in same order as cross here.
      for(i in names(newdata))
        cross$pheno[[i]] <- newdata[[i]]
    }
    else {
      ## The row.names of newdata must be index.
      mat <- match(as.character(cross$pheno[[index]]),row.names(newdata))
      if (length(mat[is.na(mat)])==length(mat)) 
        stop("no row names of newdata match index")
      cross$pheno <- cbind(cross$pheno,newdata[mat,, drop = FALSE])
    }
  }
  cross
}
