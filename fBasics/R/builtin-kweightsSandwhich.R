
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:            
#  .kweigths
################################################################################


# Code borrowed from 
#   R's contributed package "sandwich" written by Thomas Lumley
#   and Achim Zeileis.


# Rmetrics:
#   To run these functions under Rmetrics we have them implemented 
#   here as a builtin.
#   The reason is that the dependences in the original package requires 
#   zoo which may create conflicts with Rmetrics timeDate/timeSeries.


# Package: sandwich                                                       
# Version: 2.2-1                                                          
# Date: 2009-02-05                                                        
# Title: Robust Covariance Matrix Estimators                              
# Author: Thomas Lumley, Achim Zeileis                                    
# Maintainer: Achim Zeileis <Achim.Zeileis@R-project.org>                 
# Description: Model-robust standard error estimators for cross-sectional,
#   time series and longitudinal data.                         
# LazyLoad: yes                                                           
# LazyData: yes                                                           
# Depends: R (>= 2.0.0), stats, zoo                                       
# Suggests: car, lmtest, strucchange, AER, survival, MASS                 
# Imports: stats                                                          
# License: GPL-2                                                          
# Copyright: (C) 2004 Thomas Lumley and Achim Zeileis


# ------------------------------------------------------------------------------


.kweights <- 
function(x, kernel = c("Truncated", "Bartlett", "Parzen",
    "Tukey-Hanning", "Quadratic Spectral"), normalize = FALSE)
{
    kernel <- match.arg(kernel)
    if(normalize) {
      ca <- switch(kernel,  
        "Truncated" = 2,
        "Bartlett" = 2/3,
        "Parzen" = .539285,
        "Tukey-Hanning" = 3/4,
        "Quadratic Spectral" = 1)
    } else ca <- 1
    
    switch(kernel,  
        "Truncated" = { ifelse(ca * x > 1, 0, 1) },
        "Bartlett" = { ifelse(ca * x > 1, 0, 1 - abs(ca * x)) },
        "Parzen" = { 
          ifelse(ca * x > 1, 0, ifelse(ca * x < 0.5,
            1 - 6 * (ca * x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3))
        },
        "Tukey-Hanning" = {
          ifelse(ca * x > 1, 0, (1 + cos(pi * ca * x))/2)
        },
        "Quadratic Spectral" = {
          y <- 6 * pi * x/5
          ifelse(x < 1e-4, 1, 3 * (1/y)^2 * (sin(y)/y - cos(y)))
        })
}


################################################################################

