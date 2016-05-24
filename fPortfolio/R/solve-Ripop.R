
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
#  solveRipop                   Portfolio interface to solver Ripop
#  .ripopArguments              Returns arguments for solver
################################################################################


solveRipop <-
    function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Portfolio interface to solver Ripop
    
    # Example:
    #    data <- 100 * LPP2005.RET[, 1:6]
    #    spec <- portfolioSpec()
    #    setTargetReturn(spec) <- mean(data)
    #    setSolver(spec) <- "solveRipop"
    #    solveRipop(data, spec)  

    # FUNCTION:

    # Transform Data:
    Data <- portfolioData(data, spec)
    data <- getSeries(Data)
    nAssets <- getNAssets(Data)

    # Compile Arguments for Solver:
    args <- .ripopArguments(data, spec, constraints)

    # Optimize Portfolio:
    ans <- try(ipopQP(
        objective = args$objective, 
        par.lower = args$par.lower, 
        par.upper = args$par.upper, 
        eqA = args$eqA, 
        eqA.bound = args$eqA.bound, 
        ineqA = args$ineqA, 
        ineqA.lower = args$ineqA.lower, 
        ineqA.upper = args$ineqA.upper, 
        control = list()), 
        silent = TRUE)
        
    if (inherits(ans, "try-error")) {
        # When Optimization Failed:
        ans <- list(
            opt = list(dvec=NA, Dmat=NULL),
            objective = 1e99,
            status = 1,
            message = "error",
            weights = rep(0, times=nAssets))
        run <- "failed"
    } else {
        # Set Tiny Weights to Zero:
        ans$weights <- .checkWeights(ans$solution)
        ans$solution <- NULL
        run <- "passed"
    }
    ans$solver <- "solveRipop"
    ans$solution <- ans$weights
    attr(ans$weights, "invest") <- sum(ans$weights)
    attr(ans$opt, "args") <- args

    # Class:
    class(ans) <- c("solver", "list")

    # Return Value:
    ans
}
       

# -----------------------------------------------------------------------------


.ripopArguments <-
    function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns ipopQP conform arguments for the solver

    # FUNCTION:

    # Data and Constraints as S4 Objects:
    Data <- portfolioData(data)
    data <- getSeries(Data)
    nAssets <- getNAssets(Data)
    Sigma <- getSigma(Data)

    # Box Constraints:
    par.lower <- minWConstraints(data, spec, constraints)
    par.upper <- maxWConstraints(data, spec, constraints)
        
    # Set up Equality Constraints:
    eqsumW <- eqsumWConstraints(data, spec, constraints)
    eqA <- eqsumW[, -1]
    eqA.bound <- eqsumW[, 1]
    
    # Set up Inequality Constraints:
    minsumW <- minsumWConstraints(data, spec, constraints)
    maxsumW <- maxsumWConstraints(data, spec, constraints)
    ineqA <- NULL
    if(!is.null(minsumW)) ineqA = rbind(ineqA, minsumW[, -1])
    if(!is.null(maxsumW)) ineqA = rbind(ineqA, maxsumW[, -1])
    ineqA.lower <- ineqA.upper <- NULL
    if(!is.null(minsumW)) {
        ineqA.lower = c(ineqA.lower, +minsumW[, 1])
        ineqA.upper = c(ineqA.upper, rep(sum(par.upper), times=length(minsumW[, 1]))) }
    if(!is.null(maxsumW)) {
        ineqA.lower = c(ineqA.lower, rep(sum(par.lower), times=length(maxsumW[, 1])))
        ineqA.upper = c(ineqA.upper, maxsumW[, 1]) }

    # Return Value:
    list(
        objective = list(dvec = rep(0, nAssets), Dmat = Sigma), 
        par.lower = par.lower,
        par.upper = par.upper,
        eqA = eqA,
        eqA.bound = eqA.bound ,
        ineqA = ineqA,
        ineqA.lower = ineqA.lower,
        ineqA.upper = ineqA.upper
        )
}


################################################################################

