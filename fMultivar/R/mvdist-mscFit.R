
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
# FUNCTION:           DESCRIPTION:
#  dmsc                Multivariate Skew Cauchy Density function
#  pmsc                Multivariate Skew Cauchy Probability function
#  rmsc                Multivariate Skew Cauchy Random Number generator
# REQUIRES:           DESCRIPTION:
#  sn                  Contributed R-Package
#  fDISTSFIT           fBasics Package
###############################################################################


# NOTE:
#   The former multivariate skew Cauchy distribution functions have 
#   been deprecated. Use instead the functions directly from contributed
#   Package "sn".


# NOTE:
#   The former multivariate "mvFit" parameter estimation functions have
#   been deprecated. New easy to use fitting functions have been adde:
#   ms[cdt]Fit have been added.


mscFit <- 
  function(x, trace=FALSE, title=NULL, description=NULL)
{
  fit <- sn::mst.mple(
    x = rep(1, nrow(x)), y = x, start=NULL, fixed.nu=1, 
    trace=trace, penalty=NULL)
  fit$estimated <- fit$dp
  
  if (is.null(title)) 
    title <- "Skew Cauchy Parameter Estimation"
  
  if (is.null(description)) 
    description <- description()
  
  new("fDISTFIT", 
      call = match.call(), 
      model = "Skew Cauchy Distribution", 
      data = as.data.frame(x), 
      fit = fit, 
      title = title, 
      description = description)
}


###############################################################################


