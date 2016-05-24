
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
#  dmst                Multivariate Skew Student-t Density function
#  pmst                Multivariate Skew Student-t Probability function
#  rmst                Multivariate Skew Student-t Random Number generator
# REQUIRES:           DESCRIPTION:
#  sn                  Contributed R-Package
#  fDISTSFIT           fBasics Package
################################################################################


# NOTE:
#   The former multivariate skew Student-t distribution functions have 
#   been deprecated. Use instead the functions directly from contributed
#   Package "sn".


# NOTE:
#   The former multivariate "mvFit" parameter estimation functions have
#   been deprecated. New easy to use fitting functions have been adde:
#   ms[cdt]Fit have been added.


mstFit <- 
  function(x, fixed.nu=NULL, trace=FALSE, title=NULL, description=NULL)
{
  # Fit distributional Parameters:
  if (is.null(fixed.nu)) {
    fit <- sn::mst.mple(
      x = rep(1, nrow(x)), y = x, start=NULL, fixed.nu=NULL, 
      trace=trace, penalty=NULL)
    fit$estimated <- fit$dp
  } else {
    fit <- sn::mst.mple(
      x = rep(1, nrow(x)), y = x, start=NULL, fixed.nu=fixed.nu, 
      trace=trace, penalty=NULL)
    fit$estimated <- list(fit$dp, nu=fixed.nu)
  }
  
  # Add Title and Description:
  if (is.null(title)) 
    title <- "Student-t Parameter Estimation"
  
  if (is.null(description)) 
    description <- description()
  
  # Return Value:
  new("fDISTFIT", 
      call = match.call(), 
      model = "Skew Student-t Distribution", 
      data = as.data.frame(x), 
      fit = fit, 
      title = title, 
      description = description)
}


###############################################################################


