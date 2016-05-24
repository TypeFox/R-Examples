
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
#  dmsn                Multivariate Skew Normal Density function
#  pmsn                Multivariate Skew Normal Probability function
#  rmsn                Multivariate Skew Normal Random Number generator
# REQUIRES:           DESCRIPTION:
#  sn                  Contributed R-Package
#  fDISTSFIT           fBasics Package
###############################################################################


# NOTE:
#   The former multivariate skew Normal distribution functions have 
#   been deprecated. Use instead the functions directly from contributed
#   Package "sn".


# NOTE:
#   The former multivariate "mvFit" parameter estimation functions have
#   been deprecated. New easy to use fitting functions have been adde:
#   ms[cdt]Fit have been added.


msnFit <- 
  function(x, trace=FALSE, title=NULL, description=NULL)
{
  fit <- sn::msn.mle(x = rep(1, nrow(x)), y = x, start=NULL, trace=trace)
  fit$estimated <- fit$dp
  
  if (is.null(title)) 
    title <- "Skew Normal Parameter Estimation"
  
  if (is.null(description)) 
    description <- description() 
  
  new("fDISTFIT", 
      call = match.call(), 
      model = "Skew Normal Distribution", 
      data = as.data.frame(x), 
      fit = fit, 
      title = title, 
      description = description)
}


###############################################################################


.mnFit <- 
  function(x, trace = FALSE, title = NULL, description = NULL)
  {
    fit <- list()
    fit$dp <- NA
    fit$logL <- NA
    fit$aux <- NA
    fit$opt.method <- NA
    fit$estimated <- list(
      beta = colMeans(x),
      Omega = cov(x),
      alpha = rep(0, times = ncol(x)),
      nu = Inf)
    
    if (is.null(title)) 
      title <- "Normal Parameter Estimation"
    
    if (is.null(description)) 
      description <- description()
    
    new("fDISTFIT", 
        call = match.call(), 
        model = "Skew Normal Distribution", 
        data = as.data.frame(x), 
        fit = fit, 
        title = title, 
        description = description)
  }


###############################################################################


