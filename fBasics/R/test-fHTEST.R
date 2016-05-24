
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
#  'fHTEST'                  S4 Class Representation
#  show.fHTEST               S4 Print Method
################################################################################


setClass("fHTEST",
    representation(
        call = "call",
        data = "list",
        test = "list",
        title = "character",
        description = "character")
)


# ------------------------------------------------------------------------------


setMethod("show", "fHTEST",
      function(object)
{
    # A function implemented by Diethelm Wuertz

    # Source:
    #   This function copies code from base:print.htest

    # FUNCTION:

    # Unlike print the argument for show is 'object'.
    x = object

    # Title:
    cat("\nTitle:\n ", x@title, "\n", sep = "")

    # Call:
    # cat("\nCall:\n", deparse(x@call), "\n", sep = "")

    # Data Name:
    # cat("\nData Name:\n", ans@data.name, "\n", sep = "")

    # Test Results:
    test = x@test
    cat("\nTest Results:\n", sep = "")

    # Tests from tseries package:

    # Parameter:
    if (!is.null(test$parameter)) {
        parameter = test$parameter
        Names = names(parameter)
        cat("  PARAMETER:\n")
        for ( i in 1: length(Names) )
            cat(paste("    ", names(parameter[i]), ": ",
                format(round(parameter[i], 3)), "\n", sep = "") )
    }

    # Sample Estimates:
    if (!is.null(test$estimate)) {
        estimate = test$estimate
        Names = names(estimate)
        cat("  SAMPLE ESTIMATES:\n")
        for (i in 1:length(Names)) {
            cat(paste("    ", Names[i], ": ",
                round(estimate[i], digits = 4), "\n", sep = "" ) )
        }
    }

    # Statistic:
    if (!is.null(test$statistic)) {
        statistic = test$statistic
        Names = names(statistic)
        cat("  STATISTIC:\n")
        for (i in 1:length(Names)) {
            if (!is.na(statistic[i])) {
                cat(paste("    ", Names[i], ": ",
                    round(statistic[i], digits = 4), "\n", sep = "" ) )
            }
        }
    }

    # P-Value:
    if (!is.null(test$p.value)) {
        pval = test$p.value
        Names = names(pval)
        if (Names[1] == "") space = "" else space = ": "
        cat("  P VALUE:\n")
        for (i in 1:length(Names)) {
            if (!is.na(pval[i])) {
                if (class(version) != "Sversion") {
                    cat(paste("    ", Names[i], space,
                    format.pval(pval[i], digits = 4), " \n", sep = "" ) )
                } else {
                    cat(paste("    ", Names[i], space,
                    round(pval[i], digits = 4), " \n", sep = "" ) )
                }
            }
        }
    }

    # Confidence Interval:
    if (!is.null(test$conf.int)) {
        conf = test$conf.int
        # For SPlus compatibility use dimnames istead of colnames!
        colNames = dimnames(conf)[[2]]
        cat("  CONFIDENCE INTERVAL:\n")
        for (i in 1:length(colNames)) {
            cat(paste("    ", colNames[i], ": ",
                round(conf[1, i], digits = 4), ", ",
                round(conf[2, i], digits = 4), "\n", sep = "" ) )
        }
    }

    # More Specific Output Results:
    if (!is.null(test$output)) {
        cat(test$output, fill = FALSE, sep = "\n")
    }

    # Description:
    cat("\nDescription:\n ", x@description, sep = "")
    cat("\n\n")

    # Return Value:
    #   invisible()  # made visible by DW
})


################################################################################

