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
"fastTps" <- function(x, Y, m = NULL, p = NULL, theta, 
    lon.lat = FALSE, find.trA=TRUE,lambda=0, ...) {
    x <- as.matrix(x)
    d <- ncol(x)
    if (is.null(p)) {
        if (is.null(m)) {
            m <- max(c(2, ceiling(d/2 + 0.1)))
        }
        p <- (2 * m - d)
        if (p <= 0) {
            stop(" m is too small  you must have 2*m -d >0")
        }
    }
    # special arguments to send to the wendland covariance/taper function.
    # see nearest.dist for some explanation of 'method'
    cov.args <- list(k = p, Dist.args = list(method = ifelse(!lon.lat, 
        "euclidean", "greatcircle")))
    if( lambda==0){
      warning("fastTps will interpolate observations")}
    object<-mKrig(x, Y, cov.function = "wendland.cov", m = m, cov.args = cov.args, 
        theta = theta, find.trA = find.trA,lambda=lambda, ...)
    object$call<- match.call()
     class(object) <- c( "fastTps", "mKrig")
    return( object)
}

predict.fastTps <- function(object, xnew = NULL, grid.list=NULL,
                            ynew = NULL,  derivative = 0,
                            Z = NULL, drop.Z = FALSE, just.fixed = FALSE, xy=c(1,2), ...)
  {
    # the main reason to pass new args to the covariance is to increase
    # the temp space size for sparse multiplications
    # other optional arguments from mKrig are passed along in the
    # list object$args
    cov.args <- list(...)
    # predict using grid.list or as default observation locations
    if( !is.null(grid.list)){
        xnew<- make.surface.grid( grid.list)
    }
    if( is.null(xnew) ) {
        xnew <- object$x
    }   
    if (!is.null(ynew)) {
        coef.hold <- mKrig.coef(object, ynew)
        c.coef <- coef.hold$c
        d.coef <- coef.hold$d
    }
    else {
        c.coef <- object$c
        d.coef <- object$d
    }
    # fixed part of the model this a polynomial of degree m-1
    # Tmatrix <- fields.mkpoly(xnew, m=object$m)
    #
    if (derivative == 0){
        if (drop.Z | object$nZ == 0) {
            # just evaluate polynomial and not the Z covariate
            temp1 <- fields.mkpoly(xnew, m = object$m) %*% d.coef[object$ind.drift, ]
        }
        else{
            if( is.null(Z)) {
                 Z <- object$Tmatrix[, !object$ind.drift]
            }
            temp1 <- cbind(fields.mkpoly(xnew, m = object$m), Z) %*% d.coef
        }
    }
    else{
        if (!drop.Z & object$nZ > 0) {
            stop("derivative not supported with Z covariate included")
        }
        temp1 <- fields.derivative.poly(xnew, m = object$m, d.coef[object$ind.drift, 
            ])
    }
    if (just.fixed) {
        return(temp1)
    }
 
    useFORTRAN<-  (ncol(object$x)==2) & (object$args$k == 2) & (derivative==0) & (!is.null(grid.list))
  
    
    # add nonparametric part. 
    # call FORTRAN under a specific case  
    if( useFORTRAN){
        
        temp2<- multWendlandGrid(grid.list, object$knots, delta=object$args$theta, c.coef, xy=xy)
      
    }
    else{
        temp2 <- do.call(object$cov.function.name, c(object$args, 
                list(x1 = xnew, x2 = object$knots, C = c.coef, derivative = derivative), 
                cov.args))
    }  
# add two parts together
    return((temp1 + temp2))
}

multWendlandGrid <- function( grid.list,center, delta, coef, xy= c(1,2) ){
     xGrid<- grid.list[[xy[1]]]
     yGrid<- grid.list[[xy[2]]]
     mx<- length( xGrid)
     my<- length( yGrid)
# transform centers to correspond to integer spacing of grid:
# i.e. 1:nx and 1:ny
     dx<- (xGrid[mx] - xGrid[1]) / (mx-1)
     dy<- (yGrid[my] - yGrid[1]) / (my-1)
     centerScaled<- cbind( ((center[,1] - xGrid[1]) / dx) + 1,
                           ((center[,2] - yGrid[1]) / dy) + 1 )
     deltaX<- delta/dx
     deltaY<- delta/dy
     
     nc<- nrow( center)
     out<-.Fortran( "multWendlandG", PACKAGE="fields",
                 mx=as.integer(mx),
                 my=as.integer(my),
                 deltaX= as.double( deltaX),
                 deltaY= as.double( deltaY),                  
                 nc= as.integer(nc),
                 center=as.double(centerScaled),
                 coef=as.double(coef),
                 h= as.double(matrix(0,mx,my)),
                 flag=as.integer(1)
                )
     if( out$flag!= 0){
       stop("error in multWendlandG FORTRAN")}
     return( out$h)
   }

#
#"sim.fastTps.approx"<- function(fastTpsObject,...){
#                           sim.mKrig.approx( fastTpsObject,...)}
#


"sim.fastTps.approx" <- function(fastTpsObject, predictionPointsList,
                            simulationGridList=NULL, gridRefinement=5, gridExpansion=1 + 1e-07,
    M = 1, delta=NULL,  verbose = FALSE, ... ) {
    # create grid if not passed
    if( ncol( fastTpsObject$x) != 2){
      stop("Only implemented for 2 dimensions")
    }
# coerce names of grid to be x and  y    
    names(predictionPointsList) <- c( "x", "y")
    nx<- length((predictionPointsList$x))
    ny<- length((predictionPointsList$y))
    
    simulationGridList<- makeSimulationGrid2( fastTpsObject, predictionPointsList ,
                               gridRefinement, gridExpansion)
    nxSimulation<- length(simulationGridList$x)    
    nySimulation<- length(simulationGridList$y)
    sigma <- fastTpsObject$sigma.MLE
    rho <- fastTpsObject$rho.MLE
    #
    # set up various sizes of arrays
    nObs <- nrow(fastTpsObject$x)
    if (verbose) {
        cat("nObs, sigma, rho", nObs, sigma, rho, fill = TRUE)
    }
    #
    # set up object for simulating on a grid
    #
#    print( system.time(
    covarianceObject <- wendland.image.cov( 
                            setup = TRUE, grid =simulationGridList,
                            cov.args=fastTpsObject$args )
#    ))                   
     if (verbose) {
        cat( "dim of full circulant matrix ", dim(covarianceObject$wght), fill = TRUE)
    }
    # output array
    out <- matrix(NA, nx*ny, M ) 
    #
    # find conditional mean field from initial fit
    # don't multiply by sd or add mean if this is
    # a correlation model fit.
    # (these are added at the predict step).
    # from now on all predicted values are on the grid
    # represented by a matrix
    hHat<- predict.fastTps(fastTpsObject, grid.list=predictionPointsList,...)
    # empty image object to hold simulated fields
    hTrue<- c( simulationGridList, list( z= matrix(NA, nxSimulation,nySimulation)))
    ##########################################################################################
    ### begin the big loop
    ##########################################################################################
    xData<-  fastTpsObject$x
    weightsError<- fastTpsObject$weights
    for (k in 1:M) {
        # simulate full field
        if( verbose){
          cat( k, " ")}
        hTrue$z <- sqrt(rho) * sim.rf(covarianceObject)
        #
        # NOTE: fixed part of model (null space) need not be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        #   bilinear interpolation to approximate values at data locations
        #
        hData <- interp.surface(hTrue,xData)
        hPredictionGrid<- c(interp.surface.grid(hTrue, predictionPointsList)$z)
        ySynthetic <- hData + sigma * 1/sqrt(weightsError)* rnorm(nObs)
        # predict at grid using these data
        # and subtract from synthetic 'true' value

        spatialError <- c(
                          predictSurface.fastTps(fastTpsObject, 
                             grid.list = predictionPointsList,
                                 ynew  = ySynthetic, ...
                          )$z
                          ) - hPredictionGrid
        # add the error to the actual estimate  (conditional mean)
        out[ , k] <- hHat + spatialError
    }
    return( list( predictionPointsList=predictionPointsList, Ensemble=out, call=match.call()) )
}


  makeSimulationGrid2<-function( fastTpsObject, predictionPointsList,
                               gridRefinement, gridExpansion){
         nx<- length((predictionPointsList$x))
         ny<- length((predictionPointsList$y)) 
        nxSimulation<- nx*gridRefinement*gridExpansion
        nySimulation<- ny*gridRefinement*gridExpansion
        # range should include prediction grid and the data locations 
        xRange<- range(c(fastTpsObject$x[,1], predictionPointsList$x) )
        yRange<- range(c(fastTpsObject$x[,2], predictionPointsList$y) )
        midpointX<-               (xRange[2] + xRange[1])/2
        midpointY<-               (yRange[2] + yRange[1])/2
        deltaX<-    gridExpansion*(xRange[2] - xRange[1])/2 
        deltaY<-    gridExpansion*(yRange[2] - yRange[1])/2       
        return(
           list( x= seq( midpointX - deltaX, midpointX + deltaX,, nxSimulation),
                 y= seq( midpointY - deltaY, midpointY + deltaY,, nySimulation) )
                )
}
