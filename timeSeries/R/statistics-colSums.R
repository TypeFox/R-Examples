
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 COLUMN STATISTICS:
#  colStats                  Computes sample statistics by column
#  colSums                   Computes sums of all values in each column
#  colMeans                  Computes means of all values in each column
#  colSds                    Computes standardard deviation of each column
#  colVars                   Computes sample variance by column
#  colSkewness               Computes sample skewness by column
#  colKurtosis               Computes sample kurtosis by column
#  colMaxs                   Computes maximum values in each colum
#  colMins                   Computes minimum values in each colum
#  colProds                  Computes product of all values in each colum
#  colQuantiles              Computes quantiles of all values in each colum
# DEPRECATED:               NO LONGER USED:
#  colAvgs                   Computes sample mean by column
#  colStdevs                 Computes sample standard deviation by column
#  mean.timeSeries           Computes sample means by column
#  var.timeSeries            Computes sample variance by column
################################################################################


colStats <-
    function(x, FUN, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes sample statistics by column

    # FUNCTION:

    # Statistics:
    if (inherits(x, "timeSeries"))
        apply(na.omit(getDataPart(x), ...), 2, FUN, ...) #<< YC : as.matrix is slow !
    else
        apply(na.omit(as.matrix(x), ...), 2, FUN, ...)
}


# ------------------------------------------------------------------------------


# YC important because default colSums is unefficient since it retrieves
# full dimnames, i.e. rownames which is very time consuming


if (getRversion() < "2.9.0") {
    setMethod("colSums", "timeSeries",
              function(x, na.rm = FALSE, dims = 1L)
          {
              x <- getDataPart(x)
              callGeneric()
          })
} else {
    setMethod("colSums", "timeSeries",
              function(x, na.rm = FALSE, dims = 1L, ...)
          {
              x <- getDataPart(x)
              callGeneric()
          })
}


# ------------------------------------------------------------------------------


# YC important because default colSums is unefficient since it retrieves
# full dimnames, i.e. rownames which is very time consuming


if (getRversion() < "2.9.0") {
    setMethod("colMeans", "timeSeries",
              function(x, na.rm = FALSE, dims = 1L)
          {
              x <- getDataPart(x)
              callGeneric()
          })
} else {
    setMethod("colMeans", "timeSeries",
              function(x, na.rm = FALSE, dims = 1L, ...)
          {
              x <- getDataPart(x)
              callGeneric()
          })
}

# ------------------------------------------------------------------------------


colSds <- function(x, ...) { colStats(x, "sd", ...) }


colVars <- function(x, ...) { colStats(x, "var", ...) }


colSkewness <- function(x, ...) { colStats(x, "skewness", ...) }


colKurtosis <- function(x, ...) { colStats(x, "kurtosis", ...) }


colMaxs <- function(x, ...) { colStats(x, "max", ...) }


colMins <- function(x, ...) { colStats(x, "min", ...) }


colProds <- function(x, ...) { colStats(x, "prod", ...) }


# ------------------------------------------------------------------------------


colQuantiles <-
 function(x, prob = 0.05, ...)
{
    # FUNCTION:

    stopifnot(length(prob) == 1)
    colStats(x, "quantile", probs = prob, ...)
}


################################################################################
# DEPRECATED:


colAvgs <- 
  function(x, ...) 
{ 
    # FUNCTION:
    
    colMeans(x, ...) 
}


# ------------------------------------------------------------------------------


colStdevs <- 
  function(x, ...) 
{ 
    # FUNCTION:
    
    colStats(x, "sd", ...) 
}


# ------------------------------------------------------------------------------


# mean.timeSeries <- colMeans
# var.timeSeries <- colVars


################################################################################

