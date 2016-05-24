# DESP/R/DESP_SqPartCorr.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_SqPartCorr <-
function(B,n) {
  # compute squared partial correlations
  # we keep only the elements upper than a threshold. For example, if thresh=0.01, therefore SPC[i,j] equal 0 means that the proportion of the variance of the variable i not explained by the other variables, which is explained by the variable j is less than 1%.

  thresh = min(0.01, 1/sqrt(n));

  SPC = B*t(B);
  SPC = SPC*(SPC>thresh);

  return(SPC)
}
