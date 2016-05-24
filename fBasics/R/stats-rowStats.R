
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
# FUNCTION:                 DESCRIPTION:
#  rowStats                  Computes sample statistics by row
#  rowSums                   Computes sums of values in each row
#  rowMeans                  Computes means of values in each row
#  rowSds                    Computes standard deviation of each row
#  rowVars                   Computes sample variance by row
#  rowSkewness               Computes sample skewness by row
#  rowKurtosis               Computes sample kurtosis by row
#  rowMaxs                   Computes maximum values in each row
#  rowMins                   Computes minimum values in each row
#  rowProds                  Computes product of values in each row
#  rowQuantiles              Computes product of values in each row
# FUNCTION:                 NO LONGER USED:
#  rowAvgs                   Computes sample mean by row
#  rowStdevs                 Computes sample variance by row
################################################################################


rowStats <-
function(x, FUN, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes sample statistics by row

    # FUNCTION:

    # Statistics:
    apply(na.omit(as.matrix(x), ...), 1, FUN, ...)
}


# ------------------------------------------------------------------------------


# rowSums() 
#   is part of R's base package  


# ------------------------------------------------------------------------------


# rowMeans <-
#   is part of R's base package    


# ------------------------------------------------------------------------------


rowSds <- function(x, ...) { rowStats(x, "sd", ...) }
rowVars <- function(x, ...) { rowStats(x, "var", ...) }
rowSkewness <- function(x, ...) { rowStats(x, "skewness", ...) }
rowKurtosis <- function(x, ...) { rowStats(x, "kurtosis", ...) }
rowMaxs <- function(x, ...) { rowStats(x, "max", ...) }
rowMins <- function(x, ...) { rowStats(x, "min", ...) }
rowProds <- function(x, ...) { rowStats(x, "prod", ...) }


# ------------------------------------------------------------------------------


rowQuantiles <-
function(x, prob = 0.05, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:

    # Row Quantiles:
    stopifnot(length(prob) == 1)
    rowStats(x, "quantile", probs = prob, ...)
}


################################################################################


rowAvgs <- function(x, ...) rowMeans(x, ...)
rowStdevs <- rowSds


################################################################################


