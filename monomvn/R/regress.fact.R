#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## regress.fact:
##
## fit y2 ~ y1[,1:numf]  using LM (i.e., least squares)


'regress.fact' <-
  function(y1, y2, numf=1, verb=0)
{
  ## check numf
  if(length(numf) != 1 || numf < 1)
    stop("numf must be a non-negative scalar integer, you have ", numf)
  
  ## number of regressions
  numreg <- ncol(y2); if(is.null(numreg)) numreg <- 1
  if(numf > numreg) numf <- numreg
  
  ## add to progress meter
  if(verb > 0) cat(paste("using fact ", sep=""))

  ## extract the first numf factors of y1
  numcov <- ncol(y1)
  y1 <- y1[,1:numf]
  
  ## standard least-squares regression
  reglst <- lm(y2 ~ y1)
  ## cat("(", numcov-1, ")=", reglst$coef[2], " ", sep="")
  bvec <- matrix(reglst$coef, ncol=numreg)
  bvec <- rbind(bvec, matrix(rep(0, numreg*(numcov-numf)), ncol=numreg))
  res <- matrix(reglst$resid, ncol=numreg)
  actual.method <- "fact"
  ncomp <- rep(NA, numreg)
  
  return(list(method=actual.method, ncomp=ncomp, b=bvec, res=res))
}
