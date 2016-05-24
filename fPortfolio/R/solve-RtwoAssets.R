
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                    DESCRIPTION:
#  .solveRtwoAssets             Two Assets LongOnly MV Portfolio
#  .cvarSolveTwoAssets          Two Assets LongOnly CVaR Portfolio
#  .madSolveTwoAssets           Two Assets LongOnly MAD Portfolio
################################################################################


.solveRtwoAssets <-
    function(data, spec, constraints)
{
    # Description:
    #   Two Assets LongOnly MV Portfolio

    # Details:
    # ... this is only thohgt for 'unlimited' LongOnly
    # box and group constraints are discarded here.

    # FUNCTION:

    # Solver:
    # print(".mvSolveTwoAssets")

    # Convert Data and Constraints to S4 Objects:
    Data = portfolioData(data, spec)
    data <- getSeries(Data)
    Constraints = portfolioConstraints(Data, spec, constraints)

    # Stop if the Target Return is not Defined!
    targetReturn = getTargetReturn(spec)
    stopifnot(is.numeric(targetReturn))

    # Optimize Portfolio:
    nAssets = getNAssets(Data)

    # Solve the two Assets Case Analytically:
    mu <- getMu(Data)
    Sigma <- getSigma(Data)
    stopifnot(targetReturn >= min(mu))
    stopifnot(targetReturn <= max(mu))
    weights <- (targetReturn-mu[2]) / (mu[1]-mu[2])
    weights <- c(weights, 1 - weights)

    # Output List:
    ans = list(
        solver = "MVTwoAssets",
        optim = NA,
        weights = weights,
        targetReturn = targetReturn,
        targetRisk = NA,
        objective = sqrt(weights %*% Sigma %*% weights),
        status = 0,
        message = NA)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.cvarSolveTwoAssets <-
    function(data, spec, constraints)
{
    # Description:
    #   Two Assets LongOnly CVaR Portfolio

    # Details:
    # ... this is only thohgt for 'unlimited' LongOnly
    # box and group constraints are discarded here.

    # FUNCTION:

    # Solver:
    # print(".cvarSolveTwoAssets")

    # Convert Data and Constraints to S4 Objects:
    Data <- portfolioData(data, spec)
    data <- getSeries(Data)
    Constraints = portfolioConstraints(Data, spec, constraints)

    # Stop if the Target Return is not Defined!
    targetReturn = getTargetReturn(spec)
    stopifnot(is.numeric(targetReturn))
    targetAlpha = getAlpha(spec)

    # Optimize Portfolio:
    nAssets <- getNAssets(Data)

    # Solve the two Assets Case Analytically:
    mu <- getMu(Data)
    stopifnot(targetReturn >= min(mu))
    stopifnot(targetReturn <= max(mu))
    weights <- (targetReturn-mu[2]) / (mu[1]-mu[2])
    weights <- c(weights, 1 - weights)

    optim <- list(
        VaR = .varRisk(data, weights, targetAlpha),
        CVaR = -.cvarRisk(data, weights, targetAlpha),
        targetAlpha = targetAlpha)

    ans <- list(
        solver = "CVaRTwoAssets",
        optim = optim,
        weights = weights,
        targetReturn = targetReturn,
        targetRisk = NA,
        objective = optim$CVaR,
        status = 0,
        message = NA)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.madSolveTwoAssets <-
    function(data, spec, constraints)
{
    # Description:
    #   Two Assets LongOnly MAD Portfolio

    # Details:
    # ... this is only thohgt for 'unlimited' LongOnly
    # box and group constraints are discarded here.

    # FUNCTION:

    # Convert Data and Constraints to S4 Objects:
    Data = portfolioData(data, spec)
    data <- getSeries(Data)
    Constraints <- portfolioConstraints(Data, spec, constraints)

    # Stop if the Target Return is not Defined!
    targetReturn <- getTargetReturn(spec)
    stopifnot(is.numeric(targetReturn))

    # Optimize Portfolio:
    nAssets <- getNAssets(Data)

    # Solve the two Assets Case Analytically:
    mu = getMu(Data)
    stopifnot(targetReturn >= min(mu))
    stopifnot(targetReturn <= max(mu))
    weights = (targetReturn-mu[2]) / (mu[1]-mu[2])
    weights = c(weights, 1 - weights)
    targetRisk = mean( abs( (data - colMeans(data)) %*% weights ) )

    # Output List:
    ans = list(
        solver = "MADTwoAssets",
        optim = NA,
        weights = weights,
        targetReturn = targetReturn,
        targetRisk = targetRisk,
        objective = targetRisk,
        status = 0,
        message = NA)

    # Return Value:
    ans
}


################################################################################

