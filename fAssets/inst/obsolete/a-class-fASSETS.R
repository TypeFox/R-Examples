
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
# FUNCTION:             SIMULATION AND PARAMETER ESTIMATION:
#  'fASSETS'             Class representation for "fASSETS" Objects
# FUNCTION:             DESCRIPTION:
#  show.fASSETS          S4: Print method for an object of class fASSETS
#  plot.fASSETS          S3: Plot method for an object of class fASSETS
#  summary.fASSETS       S3: Summary method for an object of class fASSETS
################################################################################


setClass("fASSETS",

    # A class implemented by Diethelm Wuertz
    
    representation(
        call = "call",              # call: The matched function call
        method = "character",       # method: One of "mn", "msn", "mst"
        model = "list",             # model: A list(mu, Omega, alpha, df)
        data = "data.frame",        # Data: The data records
        fit = "list",               # fit: Results parameter estimation
        title = "character",        # title: A short title string
        description = "character")  # description: A brief description
)


# ------------------------------------------------------------------------------



setMethod("show", "fASSETS",
    function(object)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Print Method for an object of class fASSETS

    # Arguments:
    #   x - an object of class fASSETS

    # FUNCTION:

    # Title:
    cat("\nTitle:\n")
    cat(as.character(object@title), "\n")

    # Call:
    cat("\nCall:\n")
    cat(paste(deparse(object@call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    # Model Parameters:
    cat("\nModel Parameters:\n")
    print(object@model)

    # Description:
    cat("Description:\n")
    print(object@description)
    cat("\n")

    # Return Value:
    invisible(object)
})


# ------------------------------------------------------------------------------


plot.fASSETS <- 
    function(x, which = "ask", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Descriptions:
    #   Plots a fit from an assets data set or a model

    # Arguments:
    #   x - an object of class fASSETS
    #   ... - arguments to be passed

    # Notes:
    #   Library 'sn', is version  0.32-2 (2004-03-13),
    #     (C) 1998-2004 A. Azzalini, GPL
    #   For "fMV" objects have a look in "fMultivar".

    # FUNCTION:

    # Transform to a S4 object of class "fMV":
    object = new("fMV",
        call = x@call,
        method = x@method,
        model = x@model,
        data = x@data,
        fit = x@fit,
        title = x@title,
        description =
        x@description)

    # Use plot method for objects of class "fMV"
    plot(object, which = which, xlab = "Time", ylab = "Value", ...)

    # Return value:
    invisible(x)
}


# ------------------------------------------------------------------------------


summary.fASSETS =
    function(object, which = "all", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Descriptions:
    #   Summarizes a fit from an assets data set or a model

    # Print:
    print(object, ...)

    # Plot:
    plot(object, which = which, ...)

    # Return value:
    invisible(object)
}


################################################################################

