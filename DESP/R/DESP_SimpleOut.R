# DESP/R/DESP_SimpleOut.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (version 3) as published by
#  the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

DESP_SimpleOut <- 
function(ve, method='MAD', iqr.mult=1.5, mad.mult=2.5, mad.constant=1.4826){
  # simple detection of outliers, based on their Euclidean norm

  if(method=='IQR')
    {
    # based on interquantile range
    # any observation above Q3 + 1.5 * (Q3-Q1) is considered as an outlier
    out <- which(ve > quantile(ve, probs = c(3/4)) + iqr.mult * IQR(ve))
    }
  else if(method=='MAD')
    {
    # based on median absolute deviation
    # any observation above Q2 + 2.5 * MAD is considered as an outlier
    out <- which(ve > median(ve) + mad.mult * mad(ve, constant=mad.constant))
    }
  return(out)
}
