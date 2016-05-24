
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
#  solveRglpk.MAD               Portfolio interface to solver Rglpk
#  .madRglpkArguments           Returns MAD arguments for solver
# FUNCTION:                    DESCRIPTION:
#  .rglpk.MAD                   Wrapper to solver function
################################################################################


solveRglpk.MAD <-
    function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Portfolio interface to solver Rglpk

    # FUNCTION:

    # Settings:
    Data <- portfolioData(data, spec)
    data <- getSeries(Data)
    nAssets <- getNAssets(Data)
    type <- getType(spec)

    # Compile Arguments for Solver:
    args <- .madRglpkArguments(Data, spec, constraints)

    # Solve Multiassets Portfolio:
    ans <- .rglpk.MAD(
        obj = args$obj,
        mat = args$mat,
        dir = args$dir,
        rhs = args$rhs,
        types = args$types,
        max = args$max,
        bounds = args$bounds,
        verbose = args$verbose,
        nScenarios = args$nScenarios,
        nAssets = args$nAssets,
        targetReturn = args$targetReturn,
        Alpha = args$Alpha,
        Type = args$Type)
    ans$solver <- "solveRglpk.MAD"
    
    # Return Value:
    class(ans) = c("solver", "list")
    ans
}


################################################################################


.madRglpkArguments <-
    function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns glpk conform MAD arguments for the solver

    # Details:
    #       max/min:                   obj %*% x
    #       subject to:
    #                           mat %*% x  ?=  rhs
    #                                  dir  =  "?="
    #                              upper/lower bounds
    #
    #   Rglpk_solve_LP(obj, mat, dir, rhs, types = NULL, max = FALSE,
    #       bounds = NULL, verbose = FALSE)

    # FUNCTION:

    # Settings:
    Data <- portfolioData(data, spec)
    data <- getSeries(Data)
    nAssets <- ncol(data)
    nScenarios <- nrow(data)
    series <- getDataPart(data)
    series <- series - matrix(rep(colMeans(series), times=nScenarios), byrow=TRUE, ncol=nAssets)
    targetReturn <- getTargetReturn(spec)
    Type <- getType(spec)

    # Objective Function to be Maximized:
    objNames <- c(paste("e", 1:nScenarios, sep = ""), colnames(data))
    obj <- c(rep(1/nScenarios, nScenarios), rep(0, nAssets))
    names(obj) <- objNames

    # 1:
    # The negative MAD Equation Constraints:
    #   (-diag + [Returns-mu]) %*% (es, W) <= 0
    Aneg <- cbind(-diag(nScenarios), series)
    aneg <- rep(0, nrow(Aneg))
    dneg <- rep("<=", nrow(Aneg))

    # 2:
    # The positive MAD Equation Constraints:
    #   (+diag + [Returns-mu]) %*% (es, W) >= 0
    Apos <- cbind( diag(nScenarios), series)
    apos <- rep(0, nrow(Apos))
    dpos <- rep(">=", nrow(Apos))

    # 3 and 4:
    # The A_equal Equation Constraints: A_eq %*% x == a_eq
    eqsumW <- eqsumWConstraints(Data, spec, constraints)
    Aeq <- cbind(
        matrix(0, ncol=nScenarios, nrow=nrow(eqsumW)), 
        matrix(eqsumW[, -1], ncol=nAssets) )
    aeq <- eqsumW[, 1]
    deq <- rep("==", nrow(eqsumW))

    # 5:
    # The e_s > = 0 Equation Constraints:
    Aes <- cbind(diag(nScenarios), matrix(0, nrow=nScenarios, ncol=nAssets))
    aes <- rep(0, nrow(Aes))
    des <- rep(">=", nrow(Aes))

    # 6:
    # Group Constraints: A W >= a
    minsumW <- minsumWConstraints(Data, spec, constraints)
    if (is.null(minsumW)){
        Aminsum <- aminsum <- dminsum <- NULL
    } else {
        Aminsum <- cbind(
            matrix(0, nrow=nrow(minsumW), ncol=nScenarios),
            minsumW[, -1, drop=FALSE] )
        aminsum <- minsumW[, 1]
        dminsum <- rep(">=", nrow(minsumW))
    }

    # 7:
    # Group Constraints: A W <= b
    maxsumW <- maxsumWConstraints(Data, spec, constraints)
    if (is.null(maxsumW)){
        Amaxsum <- amaxsum <- dmaxsum <- NULL
    } else {
        Amaxsum <- cbind(
            matrix(0, nrow=nrow(maxsumW), ncol=nScenarios),
            maxsumW[, -1, drop=FALSE] )
        amaxsum <- maxsumW[, 1]
        dmaxsum <- rep("<=", nrow(maxsumW))
    }

    # Putting all Together:
    mat <- rbind(Aeq, Apos, Aneg, Aes, Aminsum, Amaxsum)
    rhs <-     c(aeq, apos, aneg, aes, aminsum, amaxsum)
    dir <-     c(deq, dpos, dneg, des, dminsum, dmaxsum)

    # Box Constraints: Upper and Lower Bounds as listn required ...
    minW <- minWConstraints(Data, spec, constraints)
    maxW <- maxWConstraints(Data, spec, constraints)
    nInd <- 1:(nScenarios+nAssets)
    bounds <- list(
        lower = list(ind = nInd, val = c(rep(  0, nScenarios), minW)),
        upper = list(ind = nInd, val = c(rep(Inf, nScenarios), maxW)) )

    # What variable Types, All Continuous:
    types <- NULL

    # Should I minimize or maximize ?
    max <- FALSE

    # Return Value:
    list(
        obj = obj, mat = mat, dir = dir, rhs = rhs,
        types = types, max = max, bounds = bounds, verbose = FALSE,
        nScenarios = nScenarios, nAssets = nAssets,
        targetReturn = targetReturn, Alpha = NA, Type = Type)
}


################################################################################


.rglpk.MAD <-
    function(obj, mat, dir, rhs, types, max, bounds, verbose,
    nScenarios, nAssets, targetReturn, Alpha, Type)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Rglpk MAD Solver

    # FUNCTION:

    # Solve - use Rglpk_solve_LP:
    optim <- Rglpk::Rglpk_solve_LP(
        obj = obj,
        mat = mat,
        dir = dir,
        rhs = rhs,
        types = types,
        max = max,
        bounds = bounds,
        verbose = verbose)

    # Extract Weights:
    weights <- .checkWeights(rev(rev(optim$solution)[1:nAssets]))
    attr(weights, "invest") = sum(weights)

    # Result:
    ans <- list(
        type = Type,
        solver = "Rglpk.MAD",
        optim = optim,
        weights = weights,
        solution = weights,
        targetReturn = targetReturn,
        targetRisk = -optim$optimum,
        objective = -optim$optimum,
        status = optim$status[[1]],
        message = "NA")

    # Return Value:
    ans
}


################################################################################

