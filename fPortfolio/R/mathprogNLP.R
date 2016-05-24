
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
# FUNCTION:                DESCRIPTION:
#   .solveNLP.demo          Mean-variance portfolio demo solved by NLP
################################################################################


.solveNLP.demo <- 
  function()
  {
    # Solve Mean-Variance Portfolio:
    
    # Load Dataset
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    
    # Load Swiss Pension Fund Data:
    nAssets <- 6
    data <- 100 * LPP2005REC[, 1:nAssets]     
    mu <- colMeans(data)
    Sigma <- cov(data)
    
    # Arguments: 
    start <- rep(1, nAssets)/nAssets
    objective <- function(x) { 0.5 * (x %*% Sigma %*% x)[[1]] }
    lower <- rep(0, nAssets)
    upper <- rep(1, nAssets)
    mat <- rbind(
      budget = rep(1, times=nAssets), 
      returns = colMeans(data))
    matLower <- c(
      budget = 1, 
      return = mean(data))
    matUpper <- matLower
    linCons <- list(mat, matLower, matUpper)
    control <- list()
    
    # donlp2 Solver:
    require(Rdonlp2)
    ans <- rdonlp2NLP(start, objective, lower, upper, linCons)
    ans
    
    # solnp Solver:
    # require(Rsolnp)
    ans <- rsolnpNLP(start, objective, lower, upper, linCons)
    ans
    
    # nlminb2 Solver:
    require(Rnlminb2)
    ans <- rnlminb2NLP(start, objective, lower, upper, linCons)
    ans
  }


################################################################################

