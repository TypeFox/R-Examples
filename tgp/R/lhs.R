#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
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
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## lhs:
##
## generate a Latin Hypercube Sample of size n within the rectangle
## provided.  The default "prior" for the sample is uniform, but the
## shape and mode arguments can be used to describe a beta distribution
## in each dimension.  The actual sample is generated C-side

"lhs" <-
function(n, rect, shape=NULL, mode=NULL)
{
  ## sanity checks
  if(length(n) != 1) stop(paste("length(n) should be 1, you have", length(n)))
  if(n < 0) stop(paste("n should be positive, you have", n))
  if(n == 0) return(NULL)

  ## get and check the rectangle dimensions
  if(is.null(dim(rect))) { ncol <- length(rect); d <- 1 }
  else { ncol <- ncol(rect); d <- dim(rect)[1] }
  if(ncol != 2) stop("ncol(rect) must be 2")

  ## check the shape argument should be positive and of length
  if(!is.null(shape) && length(shape) != d && all(shape > 0))
    stop(paste("For beta sampling, length(shape) should be ", d,
               ", you have ", length(shape), ", and all positive",
               sep=""))
  if(is.null(shape)) shape <- -1

  ## check the mode argument is positive and of length 1
  if(!is.null(mode) && length(mode) != d && all(mode > 0))
    stop(paste("To specify sampling modes, length(mode) should be ", d,
               ", you have ", length(mode), ", and all positive",
               sep=""))
  if(is.null(mode)) mode <- -1
  
  ## choose a random state for the C code
  state <- sample(seq(0,999), 3)

  ## run the C code
  ll <- .C("lh_sample", 
           state = as.integer(state),
           n = as.integer(n),
           d = as.integer(d),
           rect = as.double(rect), # no need to transpose
           shape = as.double(shape),
           mode = as.double(mode),
           s = double(n*d),
           PACKAGE="tgp"
           )
  
  ## just return the samples
  return(t(matrix(ll$s, nrow=d)))
}
