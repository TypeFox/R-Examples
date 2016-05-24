#******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
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


## calc.ents:
##
## wrapper used to calculate the predictive entropies
## in C

calc.ents <- function(pmat)
  {
    n <- nrow(pmat)
    if(any(pmat < 0 | pmat > 1)) stop("bad distribution p")
    
    return(.C("calc_ents_R",
              tmat = as.double(t(pmat)),
              n = as.integer(n),
              nc = as.integer(ncol(pmat)),
              ents = double(n),
              PACKAGE="plgp")$ents)
  }


## entropy:
##
## calculate the entropy of a discrete distribution
## in p

entropy <- function(p) {
  if(any(p < 0 | p > 1)) stop("bad distribution p")
  return(-sum(p*log(p)))
}


## entropy.bvsb:
##
## calculate the entropy of a discrete distribution
## in p considering only the two highest probabilities

entropy.bvsb <- function(p) {
  if(any(p < 0 | p > 1)) stop("bad distribution p")
  if(length(p) > 2) p <- sort(p, decreasing=TRUE)[1:2]
  p <- p/sum(p)
  return(-sum(p*log(p)))
}
