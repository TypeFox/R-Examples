
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 BASIC EXTENSIONS:
#  align                     aligns time series objects by approximation
#  align.default             align default method
#  atoms                     Extracts atoms from 'timeSeries' object
#  atoms.default             atoms default method
#  attach                    attach a database to the R path
#  attach.default            attach default method
#  colnames<-                colnames<- has become a generic function
#  colnames<-.default        colnames<- default method
#  cor                       cor has become a generic function
#  cor.default               cor default method
#  cov                       var has become a generic function
#  cov.default               var default method
#  log                       log has become a generic function
#  log.default               log default method
#  outlier                   outlier added generic function
#  outlier.default           outlier default method
#  rownames<-                rownames<- has become a generic function
#  rownames<-.default        rownames<- default method
#  rank                      rank has become a generic function
#  rank.default              rank default method
#  sample                    sample has become a generic function
#  sample.default            sample default method
#  sort                      sort has become a generic function
#  sort.default              sort default method
#  stdev                     stdev added generic function
#  stdev.default             stdev default method
#  termPlot                  termPlot has become a generic function
#  termPlot.default          termPlot default method
#  var                       var has become a generic function
#  var.default               var default method
#  volatility                volatility has become a generic function
#  volatility.default        volatility default method
################################################################################

## It is bad practice to create S3 generic from function especially
## when one wants to use S4 methods -> this can lead to tricky dispatch
## problems.

## .conflicts.OK = TRUE


# ------------------------------------------------------------------------------


## align <-
## function(x, y, xout, ...)
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("align")
## }


## # ------------------------------------------------------------------------------


## align.default <-
## function(x, y, xout, method = "linear", n = 50, rule = 1, f = 0,
##     ties = mean, ...)
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Align by Approximation:
##     ans = approx(x = x, y = y, xout = xout, method = method, n = n,
##         rule = rule, f = f, ties = ties, ...)

##     # Return Value:
##     ans
## }


# ------------------------------------------------------------------------------


## atoms <-
## function(x, ...)
## {
##     # A function implemented by Diethelm WUertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("atoms")
## }


## # ------------------------------------------------------------------------------


## atoms.default <-
## function(x, ...)
## {
##     # A function implemented by Diethelm WUertz

##     # FUNCTION:

##     # Return Value:
##     invisible(x)
## }


# ------------------------------------------------------------------------------


## attach <-
## function(what, pos = 2, name = deparse(substitute(what)),
##     warn.conflicts = TRUE)
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("attach")
## }


## # ------------------------------------------------------------------------------


## attach.default <- base::attach


# ------------------------------------------------------------------------------


## "colnames<-" =
## function(x, value)
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("colnames<-")
## }


## # ------------------------------------------------------------------------------


## `colnames<-.default` <- base::`colnames<-`


# ------------------------------------------------------------------------------


## cor <-
## function(x, y = NULL, use = "all.obs",
##     method = c("pearson", "kendall", "spearman"))
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("cor")
## }


# ------------------------------------------------------------------------------


## cor.default <- stats::cor


## # ------------------------------------------------------------------------------


## cov <-
## function(x, y = NULL, use = "all.obs",
##     method = c("pearson", "kendall", "spearman"))
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("cov")
## }


## # ------------------------------------------------------------------------------


## cov.default <- stats::cov


## # ------------------------------------------------------------------------------


## log <-
## function(x, base = exp(1))
## {   # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("log")
## }


## # ------------------------------------------------------------------------------


## log.default <-
## function(x, base = exp(1))
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     base::log(x, base)
## }


# ------------------------------------------------------------------------------


## rank <-
## function(x, na.last = TRUE,
##     ties.method = c("average", "first", "random", "max", "min"))
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("rank")
## }


# ------------------------------------------------------------------------------


## rank.default <-
## function(x, na.last = TRUE,
##     ties.method = c("average", "first", "random", "max", "min"))
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     base::rank(x, na.last = na.last, ties.method = ties.method)
## }


# ------------------------------------------------------------------------------


## sample <-
## function(x, ...)
## {   # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("sample")
## }


## # ------------------------------------------------------------------------------


## sample.default <-
## function(x, size, replace = FALSE, prob = NULL, ...)
## {
##     # FUNCTION:

##     base::sample(x, size, replace = replace, prob = prob)
## }


# ------------------------------------------------------------------------------


## if(getRversion() < "2.4.0") {

##     # Note:
##     # sort() has been S3 generic in 'base' since 2.4.0
##     # Otherwise use something that works here

##     sort <- function(x, decreasing = FALSE, ...)
##     {
##         if (!is.logical(decreasing) || length(decreasing) != 1)
##             stop("'decreasing' must be a length-1 logical vector.\nDid you intend to set 'partial'?")
##         UseMethod("sort")
##     }

##     sort.default <- function(x, decreasing = FALSE, ...) {
##         if (is.object(x))
##         x[order(x, decreasing = decreasing)]
##         else base::sort(x, decreasing = decreasing, ...)
##     }

## }# endif {only for outdated R}


# ------------------------------------------------------------------------------


## outlier <-
## function(x, sd = 5, complement = TRUE, ...)
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     UseMethod("outlier")
## }


# ------------------------------------------------------------------------------


## outlier.default <-
## function(x, sd = 5, complement = TRUE, ...)
## {
##     # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Returns outliers

##     # Arguments:
##     #   x - a numeric vector
##     #   sd - a numeric value of standard deviations, e.g. 5
##     #       means that values larger or smaller tahn five
##     #       times the standard deviation of the series will
##     #       be detected.
##     #   complement - a logical flag, should the outlier series
##     #       or its complements be returned.

##     # Note:
##     #   This function is thought to find splits in financial
##     #   price or index data. If a price or index is splitted we
##     #   observe in the returns a big jump of several standard
##     #   deviations.

##     # FUNCTION:

##     # Find Outliers:
##     SD = sd * sd(x)
##     if (complement) {
##         ans  = x[x <= SD]
##     } else {
##         ans = x[x > SD]
##         names(ans) = as.character(which(x > SD))
##     }

##     # Return Value:
##     ans
## }


# ------------------------------------------------------------------------------


## "rownames<-" =
## function(x, value)
## {   # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("rownames<-")
## }


## # ------------------------------------------------------------------------------


## `rownames<-.default` <- base::`rownames<-`


# ------------------------------------------------------------------------------


stdev.default <-
function(x, na.rm = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    stats::sd(x = x, na.rm = na.rm)
}


# ------------------------------------------------------------------------------


stdev <- function(x, na.rm = FALSE)
    UseMethod("stdev")

# ------------------------------------------------------------------------------


termPlot <-
function(model, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    UseMethod("termPlot")
}


# ------------------------------------------------------------------------------


termPlot.default <-
function(model, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    stats::termplot(model, ...)
}


# ------------------------------------------------------------------------------


## var <-
## function(x, y = NULL, na.rm = FALSE, use)
## {
##     # A function implemented by Diethelm Wuertz

##     # FUNCTION:

##     # Return Value:
##     UseMethod("var")
## }


## # ------------------------------------------------------------------------------


## var.default <- stats::var


# ------------------------------------------------------------------------------

volatility <-
function(object, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    UseMethod("volatility")
}


# ------------------------------------------------------------------------------


volatility.default <-
function(object, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Compute volatility:
    x = object
    ans = (x - mean(x))^2

    # Return Value:
    ans
}


################################################################################

