
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


###############################################################################
# FUMCTION:             DESCRIPTION REGRESSION METHODS:
#  fitted.fREG           Fitted values method for an object of class fREG
###############################################################################


setMethod(f = "fitted", signature(object = "fREG"), definition =
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fitted values method for an object of class fREG

    # FUNCTION:

    # Fitted Values:
    fitted <- object@fitted

    # Get original time series class:
    data = slot(object, "data")$data
    dataClass = class(data)[1]

    # Transform:
    if (dataClass == "timeSeries") {
        ans <- data
        data.mat <- matrix(fitted)
        rownames(data.mat) <- rownames(data)
        colnames(data.mat) <- object@data$unit
        series(ans) <- data.mat
        colnames(ans) <- as.character(object@formula[2])
    } else {
        ans <- data
    }

    # Return Value:
    ans
})


###############################################################################


