
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
# FUNCTION:                     DESCRIPTION:
#  solveRquadprog.CLA            Portfolio interface to solver Rquadprog
#  .claRquadprogArguments        Returns arguments for solver
# FUNCTION:                     DESCRIPTION:
#  .quadprog.CLA                 Wrapper to solver function
################################################################################


solveRquadprog.CLA <-
    function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Portfolio interface to solver Rquadprog

    # Example:
    #   solveRquadprogCLA(data, spec, constraints)[-3]
   
    # FUNCTION:   

    # Update Specification:
    setTargetReturn(spec) <- NA
    
    # Transform Data:
    Data <- portfolioData(data, spec)
    nAssets <- getNAssets(Data)
        
    # Compile Arguments for Solver:
    args <- .claRquadprogArguments(data, spec, constraints)
    
    # Solve Multiassets Portfolio:
    ans <- .quadprog.CLA(
        Dmat = args$Dmat, 
        dvec = args$dvec, 
        Amat = args$Amat, 
        bvec = args$bvec, 
        meq = args$meq,
        lambda = args$lambda)
     
    # Save Arguments:
    ans$optim$args <- args
    class(ans) <- c("solver", "list")

    # Return Value:
    ans
}


################################################################################


.claRquadprogArguments <-
    function(data, spec, constraints)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns quadprog conform arguments for the solver
    
    # FUNCTION:
    
    # Set up the default quadprog QP
    ans <- .rquadprogArguments(data, spec, constraints)
    
    # Optimize:
    # min(-d^T x + 1/2 x^T D x) 
    
    # Start from it and modify for CLA:
    lambda <- spec@model$param$lambda
    ans$Dmat <- lambda * getSigma(portfolioData(data)) / 2
    ans$dvec <- getMu(portfolioData(data)) 
    ans$Amat <- ans$Amat[, -1]
    ans$bvec <- ans$bvec[-1]
    ans$meq <- ans$meq - 1
    ans$dir <- ans$dir[-1]
    ans$lambda <- lambda

    # Return Value:
    ans
}


################################################################################


.quadprog.CLA <-
    function(Dmat, dvec, Amat, bvec, meq, lambda)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Goldfarb and Idnani's quadprog solver function
    
    # Note:
    #   Requires to load contributed R package quadprog from which we use
    #   the Fortran subroutine of the quadratic solver.
    
    # Package: quadprog
    #   Title: Functions to solve Quadratic Programming Problems.
    #   Author: S original by Berwin A. Turlach <berwin.turlach@anu.edu.au>
    #       R port by Andreas Weingessel <Andreas.Weingessel@ci.tuwien.ac.at>
    #   Maintainer: Andreas Weingessel <Andreas.Weingessel@ci.tuwien.ac.at>
    #   Description: This package contains routines and documentation for
    #       solving quadratic programming problems.
    #   License: GPL-2
    
    # Value of slove.QP():
    #   solution - vector containing the solution of the quadratic
    #       programming problem.
    #   value - scalar, the value of the quadratic function at the
    #       solution
    #   unconstrained.solution - vector containing the unconstrained
    #       minimizer of the quadratic function.
    #   iterations - vector of length 2, the first component contains
    #       the number of iterations the algorithm needed, the second
    #       indicates how often constraints became inactive after
    #       becoming active first. vector with the indices of the
    #       active constraints at the solution.

    # FUNCION:

    # Optimize:
    optim <- quadprog::solve.QP(
        Dmat = Dmat, 
        dvec = dvec, 
        Amat = Amat, 
        bvec = bvec, 
        meq = meq, 
        factorized = FALSE) 
    
    # Set Tiny Weights to Zero:
    weights <- .checkWeights(optim$solution)
    attr(weights, "invest") <- sum(weights) 
    
    # Compose Output List:
    Sigma <- Dmat * 2 / lambda
    mu <- dvec
    ans <- list(
        type = "MV",
        solver = "solveRquadprog.CLA",
        optim = optim,
        weights = weights,
        solution = weights,
        targetReturn = (optim$solution %*% mu)[[1,1]],
        targetRisk = sqrt(optim$solution %*% Sigma %*% optim$solution)[[1,1]],
        objective = optim$crval,
        # To do: Add status information
        status = 0,
        message = "minRisk")
            
    # Return Value:
    ans
}


################################################################################

