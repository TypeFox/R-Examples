#
# file: rmvnorm90ci_exact.R
#
# This file is part of the R-package decisionSupport
# 
# Authors: 
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF) 
#	http://www.worldagroforestry.org
# 
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
##############################################################################################
# rmvnorm90ci_exact( n, lower, upper, correlationMatrix)
##############################################################################################
#' 90\%-confidence interval multivariate normal random number generation. 
#' 
#' This function generates normally distributed multivariate random numbers which parameters are 
#' determined by the 90\%-confidence interval. The calculation of \code{mean} and \code{sd} is 
#' exact.
#' @param n \code{integer}: Number of observations to be generated.
#' @param lower \code{numeric vector}: lower bound of the 90\% confidence interval.
#' @param upper \code{numeric vector}: upper bound of the 90\% confidence interval.
#' @param correlationMatrix \code{numeric matrix}: symmetric matrix which is the correlation matrix of the 
#'   multivariate normal distribution. In particular, all diagonal elements must be equal to 1. 
#' @seealso \code{\link{random}}, \code{\link[mvtnorm]{rmvnorm}}
#' @export
rmvnorm90ci_exact <- function(n, lower, upper, correlationMatrix){
  correlationMatrix<-as.matrix(correlationMatrix)
  lower<-data.matrix(lower)
  upper<-data.matrix(upper)
  colnames(lower)<-NULL
  colnames(upper)<-NULL
  # Constants:
  # 95%-critical value of standard normal distribution (c_0.95=1.645):
  c_0.95=qnorm(0.95)
  # Check preconditions
  if ( !is.numeric(lower) || !is.numeric(upper) )
    stop("lower and upper value of the 90%-confidence interval must be given.")
  if( !identical(length(lower), length(upper)) )
    stop("lower and upper vectors must be of the same length.")
  if( !identical( correlationMatrix, t(correlationMatrix) ) ) 
    stop("correlationMatrix must be a symmetric matrix.")
  if( !identical( as.vector(diag(correlationMatrix)), rep(1, nrow(correlationMatrix)) ) )
  		stop("All diagonal elements of correlationMatrix must be equal to 1.")
  if( !identical(length(lower), nrow(correlationMatrix)) )
    stop("confidence interval vectors and correlationMatrix must the same number of rows.")
  # Calculate the mean vector from the 90-% confidence interval:
  mean<-rowMeans(cbind(lower,upper))
  # Calculate the standard deviation from the 90-% confidence interval: 
  sd<-((mean - lower)/c_0.95)[,1]
  # Calculate the covariance matrix:
  sigma<-(t(sd*correlationMatrix))*sd
  # Generate the random numbers 
  x<-mvtnorm::rmvnorm(n=n,
             mean=mean,
             sigma=sigma)
  # Return the generated random numbers:
  x
}
