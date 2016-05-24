# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"sim.mKrig.approx" <- function(mKrigObject, predictionPoints = NULL, 
    predictionPointsList = NULL, simulationGridList = NULL, gridRefinement = 5, 
    gridExpansion = 1 + 1e-07, M = 1, nx = 40, ny = 40, nxSimulation = NULL, 
    nySimulation = NULL, delta = NULL, verbose = FALSE,...) {
    if (ncol(mKrigObject$x) != 2) {
        stop("conditional simulation only implemented for 2 dimensions")
    }
    # create prediction set of points based on what is passed
    if (is.null(predictionPoints)) {
        predictionPoints <- makePredictionPoints(mKrigObject, 
            nx, ny, predictionPointsList)
    }
    if (is.null(simulationGridList)) {
        simulationGridList <- makeSimulationGrid(mKrigObject, 
            predictionPoints, nx, ny, nxSimulation, nySimulation, 
            gridRefinement, gridExpansion)
    }
    nxSimulation <- length(simulationGridList$x)
    nySimulation <- length(simulationGridList$y)
    sigma <- mKrigObject$sigma.MLE
    rho <- mKrigObject$rho.MLE
    #
    # set up various sizes of arrays
    nObs <- nrow(mKrigObject$x)
    if (verbose) {       
        cat("nObs, sigma, rho", nObs, sigma, rho, fill = TRUE)
        cat("simulationGridList)", fill=TRUE)
        print( t( stats( simulationGridList)))
    }
    # set up object for simulating on a grid using circulant embedding
    covarianceObject <- stationary.image.cov(setup = TRUE, grid = simulationGridList, 
        cov.function = mKrigObject$cov.function, cov.args = mKrigObject$args, 
        delta = delta)
    if (verbose) {
        cat("dim of full circulant matrix ", dim(covarianceObject$wght), 
            fill = TRUE)
    }
    #
    # find conditional mean field from initial fit
    hHat <- predict(mKrigObject, xnew = predictionPoints, grid.list = predictionPointsList, ...)
    # setup output array to hold ensemble
    out <- matrix(NA, length(hHat), M)
    # empty image object to hold simulated fields
    hTrue <- c(simulationGridList, list(z = matrix(NA, nxSimulation, 
        nySimulation)))
    ##########################################################################################
    ### begin the big loop
    ##########################################################################################
    xData <- mKrigObject$x
    weightsError <- mKrigObject$weights
    for (k in 1:M) {
        # simulate full field
        if (verbose) {
            cat(k, " ")
        }
        hTrue$z <- sqrt(rho) * sim.rf(covarianceObject)
        #
        # NOTE: fixed part of model (null space) need not be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        #   bilinear interpolation to approximate values at data locations
        #
        hData <- interp.surface(hTrue, xData)
        hPredictionGrid <- interp.surface(hTrue, predictionPoints)
        ySynthetic <- hData + sigma * 1/sqrt(weightsError) * 
            rnorm(nObs)
        if (verbose) {
            cat("stats for synthetic values", fill = TRUE)
            print(t(stats(ySynthetic)))
        }
        # predict at grid using these data
        # and subtract from synthetic 'true' value
        spatialError <- predict(mKrigObject, xnew = predictionPoints, 
            grid.list = predictionPointsList, ynew = ySynthetic, ...) - 
            hPredictionGrid
        # add the error to the actual estimate  (conditional mean)
        out[, k] <- hHat + spatialError
    }
    return(list(predictionPoints = predictionPoints, Ensemble = out, 
        call = match.call()))
}
makeSimulationGrid <- function(mKrigObject, predictionPoints, 
    nx, ny, nxSimulation, nySimulation, gridRefinement, gridExpansion) {
    # if prediction grid is passed use these to deterimine the simulation grid.
    #
    if (is.null(nxSimulation) | is.null(nySimulation)) {
        nxSimulation <- nx * gridRefinement * gridExpansion
        nySimulation <- ny * gridRefinement * gridExpansion
    }
    # Note NULL values are transparent ther because of 'c' operator.
    xRange <- range(c(mKrigObject$x[, 1], predictionPoints[, 
        1]))
    yRange <- range(c(mKrigObject$x[, 2], predictionPoints[, 
        2]))
    midpointX <- (xRange[2] + xRange[1])/2
    midpointY <- (yRange[2] + yRange[1])/2
    deltaX <- gridExpansion * (xRange[2] - xRange[1])/2
    deltaY <- gridExpansion * (yRange[2] - yRange[1])/2
    return(list(x = seq(midpointX - deltaX, midpointX + deltaX, 
        , nxSimulation), y = seq(midpointY - deltaY, midpointY + 
        deltaY, , nySimulation)))
}
makePredictionPoints <- function(mKrigObject, nx, 
    ny, predictionPointsList) {
    if (is.null(predictionPointsList)) {
        predictionPointsList <- fields.x.to.grid(mKrigObject$x, 
            nx = nx, ny = ny)
    }
    return(make.surface.grid(predictionPointsList))
}
