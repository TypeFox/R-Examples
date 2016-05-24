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


## addy:
##
## Use ML to estimate the mean and variance matrix of a
## random vector y=y1,y2 when there is a monotone pattern of
## missing data.  Assume that we have the mean of y1 in m1
## and the variance of y1 in s11, with the variance-covariance
## of the mean vector in c11.  We use the complete cases to
## compute the regression of y2 on y1, then update the mean
## and variance matrix accordingly.  Apply this recursively
## and you can estimate the entire mean vector (with its
## variance matrix) and variance matrix (without its variance
## matrix, unfortunately).
##
## adapted from Daniel F. Heitjan, 03.02.13


`addy` <-
function(y1, y2, m1, s11, method="plsr", p=1.0, ncomp.max=Inf, validation="CV", 
         verb=0, quiet=TRUE)
  {
    ## decide what kind of regression to do and return coeffs & mean-sq resids
    reg <- regress(y1, y2, method, p, ncomp.max, validation, verb, quiet)
    
    ## separate out the intercept term from the regression coeffs
    b0 <- reg$b[1,]
    b1 <- matrix(reg$b[-1,], ncol=ncol(reg$b))

    ## print bhat and s2-hat to the screen
    if(verb >= 2) {
      cat("\n\nbhat =\n"); print(reg$b)
      cat("\ns2hat =\n"); print(reg$S)
    }
    
    ## Update the parameters

    ## mean
    if(length(m1) == 1) {
      m2 <- b0 + b1 * m1
      s21 <- matrix(b1 * s11, nrow=ncol(reg$b))
    } else {
      m2 <- b0 + t(b1) %*% m1
      s21 <- t(b1) %*% s11
    }

    ## don't actually need to invert s11 here
    ## s22 <- reg$S + s21 %*% solve(s11,t(s21))
    ## s22 <- reg$S + t(b1) %*% s11 %*% b1
    s22 <- reg$S + s21 %*% b1

    ## perhaps print new components as we go along
    if(verb >= 2) {
      cat("\nnew components of mu: "); cat(paste(round(m2,2))); cat("\n")
      cat("\nnew cols of S:\n")
      print(rbind(t(s21), s22))
      if(verb >= 3) readline("\npress RETURN to continue: ")      
    }
    
    ## return
    return(list(method=rep(reg$method, ncol(reg$b)), ncomp=reg$ncomp,
                lambda=reg$lambda, mu=m2, s21=s21, s22=s22))
  }

