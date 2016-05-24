### fisher.g.test.R (2008-06-03)
###
###     Fisher's exact g test
###
### Copyright 2003-08 Konstantinos Fokianos and Korbinian Strimmer
###
###
### This file is part of the `GeneCycle' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


# The following function calculates the p-value of
# Fisher's exact g test given a single time series
#
# note that constant times series result in a p-value of 1
# (i.e. the null hypothesis of a purely random process is not rejected)

# Fishers exact g test (single time series, x is a vector)
fisher.g.test.single <- function(x, ...) 
{
    # constant time series result in a p-value of 1
    if( is.constant.single(x) ) return(1)
    
    f.spec <- periodogram.spec(x, ...)

    # we have to exclude the intensity at frequency \pi for 
    # an even number of time points 
    # (the need for this was pointed out by Sue Wilson)
    # changed 3 June 2008
    
    if ( (length(x) %% 2) == 0 ) 
      f.spec = f.spec[-length(f.spec)]  # remove peak at frequency \pi

    # Max Periodogram at Frequency w1 in radians/unit time:
    w1 <- which.max(f.spec)
    fisher   <- f.spec[w1]/sum(f.spec)
    upper    <- floor(1/fisher)
    compose  <- rep(NA, length=upper)
    
    m = length(f.spec)

    for (j in 1:upper)
    {
      # original code
      #compose[j]  <- (gamma(m+1)/(gamma(j+1)*gamma(m-j+1)))*((-1)^(j-1))*(1-j*fisher)^(m-1)
      
      # problem: the gamma function diverges for times series of length > 341
      # this bug and the following solution (which allows to compute g-statistic 
      # for long time series) was kindly suggested by Benjamin Tyner
      
      compose[j] <- (-1)^(j-1)*exp(lchoose(m,j)+(m-1)*log(1-j*fisher))

    }
    pval  <- sum(compose)  
    if (pval > 1) pval <- 1 # this may happen due to numerical error
    
    return(pval)
}

# Fishers exact g test (multiple time series)
fisher.g.test <- function(x, ...) 
{
  xm <- as.matrix(x)
  
  num.series <- dim(xm)[2] # number of columns
  pvalues <- rep(NA, length=num.series)
  for (i in 1:num.series)
  {
     pvalues[i] <- fisher.g.test.single(xm[,i], ...)
  }
  
  return(pvalues)
}
