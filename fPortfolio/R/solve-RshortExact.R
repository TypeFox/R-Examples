
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
#  solveRshortExact             Portfolio interface to solver RshortExact
#  .rshortExactArguments        Returns arguments for solver
#  .rshortExact                 Wrapper to solver function
################################################################################


solveRshortExact <-
    function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Portfolio interface to solver RshortExact

    # Details:
    #   If getTargetReturn() is specified we minimze the risk,
    #   if getTargetRisk() is pecified we maximize the risk.

    # FUNCTION:

    # Convert Data and Constraints to S4 Objects:
    Data <- portfolioData(data, spec)
    data <- getSeries(Data)
    Sigma <- getSigma(Data)
    Constraints <- portfolioConstraints(Data, spec, constraints)

    # Stop if the Target Return is not Defined!
    optimize <- getOptimize(spec)
    targetReturn <- getTargetReturn(spec)
    targetRisk <- getTargetRisk(spec)

    # Get '.rshortexact' conform arguments:
    args <- .rshortExactArguments(Data, spec, Constraints)

    # Solve Portfolio:
    ans <- .rshortExact(
        optimize = optimize,
        C0 = args$C0, 
        a = args$a, 
        b = args$b, 
        c = args$c, 
        d = args$d,
        Sigma = Sigma,
        invSigma = args$invSigma, 
        mu = args$mu,
        targetReturn, 
        targetRisk)

    # Return Value:
    class(ans) <- c("solver", "list")
    ans
}


# ------------------------------------------------------------------------------


.rshortExactArguments <-
function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns 'shortexact' conform arguments for the solver

    # FUNCTION:

    # Data as S4 Objects:
    Data <- portfolioData(data, spec)
    data <- getSeries(Data)

    # Get Specifications:
    mu <- getMu(Data)
    Sigma <- getSigma(Data)
    weights <- getWeights(spec)
    targetReturn <- getTargetReturn(spec)
    targetRisk <- getTargetRisk(spec)

    # Parameter Settings:
    C0 <- 1
    one <- rep(1, times = length(mu))
    invSigma <- solve(Sigma)
    a <- as.numeric(mu %*% invSigma %*% mu)
    b <- as.numeric(mu %*% invSigma %*% one)
    c <- as.numeric(one %*% invSigma %*% one)
    d <- as.numeric(a*c - b^2)

    # Return Value:
    list(C0 = C0, a = a, b = b, c = c, d = d, mu = mu, invSigma = invSigma)
}


################################################################################


.rshortExact <-
    function(optimize, C0, a, b, c, d, Sigma, invSigma, mu, 
    targetReturn, targetRisk)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Analytical 'shortexact' solver function

    # FUNCTION:

    # Optimize:
    if (!is.null(targetReturn)) {
        if (optimize == "minRisk") {
            # Compute Target Risk:
            objective <- targetRisk <-
                sqrt((c*targetReturn^2 - 2*b*C0*targetReturn + a*C0^2) / d)
        } else if (optimize == "maxReturn")  {
            # Compute Target Return:
            aq <- c
            bq <- -2*b*C0
            cq <- a*C0^2 - d*targetRisk^2
            objective <- targetReturn <- (-bq + sqrt(bq^2 - 4*aq*cq)) / (2*aq)
        }
    }
    
    # Compute Weights:
    if (is.null(targetReturn)) {
        # global minimum Variance Portfolio:
        optimize <- "minvariancePortfolio"
        one <- rep(1, times=length(mu))
        weights <- as.vector(invSigma %*% one) / ( one %*% invSigma %*% one )
        objective <- t(weights) %*% Sigma %*% weights
    } else {
        weights <- as.vector(invSigma %*% ((a-b*mu)*C0 + (c*mu-b)*targetReturn )/d)
    }
    weights <- .checkWeights(weights)

    # Return Value:
    list(
        type = "MV",
        solver = "solveRshortExact",
        optim = NA,
        solution = weights,
        weights = weights,
        targetReturn = targetReturn,
        targetRisk = targetRisk,
        objective = objective,
        status = 0,
        message = optimize)
}


################################################################################

