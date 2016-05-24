
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

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               DESCRIPTION:
#  .modelSeries            Models a timeSeries object to use formulas
################################################################################


.modelSeries <-
function(formula, data, fake = FALSE, lhs = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   formula -
    #   data - a timeSeries, a data.frame or a numeric vector
    #   fake -
    #   lhs -

    # Details:
    #   Time Series Modelling
    #   Regression Modelling
    #   Data Management

    # FUNCTION:

    # If no respnonse is pecified:
    if (length(formula) == 2) {
        formula = as.formula(paste("x", formula[1], formula[2], collapse = ""))
        stopifnot(!missing(data))
    }

    # Isf data is missing, take the first data set from the search path:
    if (missing(data)) {
        data = eval(parse(text = search()[2]), parent.frame())
    }

    if (is.numeric(data)) {
        data = data.frame(data)
        colnames(data) = all.vars(formula)[1]
        lhs = TRUE
    }

    # If we consider a faked formula:
    if (fake) {
        response = as.character(formula)[2]
        Call = as.character(match.call()[[2]][[3]])
        method = Call[1]
        predictors = Call[2]
        formula = as.formula(paste(response, "~", predictors))
    }

    # If only left hand side is required:
    if (lhs) {
        response = as.character(formula)[2]
        formula = as.formula(paste(response, "~", 1))
    }

    # Create Model Data:
    x = model.frame(formula, data)

    # Convert:
    if (class(data) == "timeSeries") x = timeSeries(x)
    if (fake) attr(x, "control") <- method

    # Return value:
    x

}


################################################################################
