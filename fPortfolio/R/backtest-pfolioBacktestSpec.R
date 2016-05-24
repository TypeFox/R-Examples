
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
# FUNCTION:                    DESCRIPTION:
#  portfolioBacktest            Returns an object of class fPFOLIOBACKTEST
################################################################################


portfolioBacktest <-
function(
    windows = list(
        windows = "equidistWindows",
        params = list(
            horizon = "12m")),
    strategy = list(
        strategy = "tangencyStrategy",
        params = list()),
    smoother = list(
        smoother = "emaSmoother",
        params = list(
            doubleSmoothing = TRUE,
            lambda = "3m",
            skip = 0,
            initialWeights = NULL)),
    messages = list() )
{
    # A function implemented by William Chen and Diethelm Wuertz
        
    # Description:
    #   Specifies a portfolio to be optimized from scratch

    # Example:
    #   portfolioBacktest()
   
    # Arguments:
    #   windows - rolling windows slot:
    #       windows - the name of the rollings windows function
    #       params - parameter list for windows settings:
    #           horizon - length of the rolling windows
    #   strategy - portfolio strategy slot:
    #       strategy - the name of the portfolio strategy function
    #       params - parameter list for strategy settings:
    #   smoother - smoother approach slot:
    #       smoother - the name of the portfolio weights smoother function
    #       params - parameter list for smoother settings:
    #           doubleSmoothing - a flag sould we double smooth the weights?
    #           lambda - length of the ema smoothing parameter
    #           skip - hoqw many periods should be skipped for smoothing ?
    #           initialWeights - vector of initial weights
   
    # FUNCTION:

    # Return Value:
    new("fPFOLIOBACKTEST",
        windows = windows,
        strategy = strategy,
        smoother = smoother,
        messages = messages)
}


################################################################################

