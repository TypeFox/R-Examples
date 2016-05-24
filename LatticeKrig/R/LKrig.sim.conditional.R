# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

LKrig.sim.conditional <- function(LKrigObj, M = 1, x.grid = NULL, 
    grid.list = NULL, nx = 80, ny = 80, ..., Z.grid = NULL, seed=42, verbose=FALSE) {
    # generate grid if not specified
    if (is.null(x.grid)) {
        if (is.null(grid.list)) {
            grid.list <- fields.x.to.grid(LKrigObj$x, nx = nx, ny = ny)
        }
        x.grid <- make.surface.grid(grid.list)
    }
    if( verbose){
    	cat("LKrig.sim.conditional: x.grid")
    	print( x.grid)
    }
    # NOTE: the name x.grid may be misleading because it just needs to a column matrix of
    # locations. It need not follow any regular pattern.
    # now generate the error surfaces
    # begin block
    # create vector of seeds if needed
    if( length(seed)==1){
        seeds<- seed + ( 0:(M-1))
      }
    #
    g.conditional.draw <-    matrix(NA, ncol = M, nrow = nrow(x.grid))

    N <- nrow(LKrigObj$y)  
    # complete set of locations to evaluate the field must
    # include the observations too
    PHIGrid<- LKrig.basis(x.grid,LKrigObj$LKinfo) 
    if( verbose){
    	cat("LKrig.sim.conditional: dim(PHIGrid)",  dim(PHIGrid), fill=TRUE)
    }
   
    # predicted field at grid from the actual data
    spatialPart<- (PHIGrid%*% LKrigObj$c.coef)
    ghat <-  spatialPart
    if( !is.null(LKrigObj$LKinfo$fixedFunction) ){
      d.coef.draw<- matrix(NA, ncol = M, nrow = length( LKrigObj$d.coef) )
   }
   else{
   	d.coef.draw<- matrix(NA, ncol = M, nrow=1)
   	}
   	ghat<- predict( LKrigObj, x= x.grid, Z = Z.grid )
    for (k in 1:M) {
        cat(k, " ")
        out<- simConditionalDraw( k, LKrigObj, ghat, x.grid, Z.grid,  PHIGrid, seeds, ..., verbose=verbose)
        if( !is.null(LKrigObj$LKinfo$fixedFunction) ){
        d.coef.draw[,k] <- out$d.coef
        }
        g.conditional.draw[, k] <- out$g.conditional
    }
    #
    return(list(x.grid = x.grid, ghat = ghat, g.draw = g.conditional.draw,
                           d.coef.draw= d.coef.draw))
}

simConditionalDraw <- function(index=1,  LKrigObj, ghat, x.grid, Z.grid, PHIGrid, seeds= 123,  verbose=FALSE){
require(LatticeKrig)
        set.seed( seeds[index] )
# generate process at grid and also on the observation locations.
        simCoefficients<- LKrig.sim(LKinfo = LKrigObj$LKinfo, just.coefficients=TRUE)
        g.unconditional.data <-LKrigObj$wX %*%simCoefficients 
        g.unconditional.data <- sqrt(LKrigObj$rho.MLE) * g.unconditional.data
        g.unconditional.grid <-sqrt(LKrigObj$rho.MLE) *PHIGrid%*%simCoefficients 
        # generate a synthetic data set with fixed part set to zero.
        N<- length( LKrigObj$y)
        y.synthetic.data <- g.unconditional.data + LKrigObj$sigma.MLE * 
            rnorm(N)
# this may confusing. divide by  sqrt(weights) to cancel out this term in the
# wX matrix for  g.unconditional.data and to adjust measurement error variance         
        y.synthetic.data<- y.synthetic.data / sqrt(LKrigObj$weights)
        # use LKrig to find the predictions for the xgrid locations
        # NOTE that LKrig will still estimate the fixed part.
        # and it is important to include this part of estimate
        obj.fit.synthetic <- LKrig(LKrigObj$x, y.synthetic.data,
                                   LKinfo = LKrigObj$LKinfo,
                                       wX = LKrigObj$wX,
                                       wU = LKrigObj$wU,
                                   lambda = LKrigObj$lambda,
                                        Z = LKrigObj$Z,
                                  weights = LKrigObj$weights,
                             use.cholesky = LKrigObj$Mc)      
        #
        # predict field
        spatialPart<- (PHIGrid%*% obj.fit.synthetic$c.coef)
        ghat.synthetic<-  spatialPart
        if( !is.null(LKrigObj$LKinfo$fixedFunction) ){       	
                 fixedPart<- predict(
                 obj.fit.synthetic, xnew=x.grid, Znew = Z.grid, just.fixed=TRUE)
                 d.coef <- obj.fit.synthetic$d.coef
                 ghat.synthetic<- ghat.synthetic + fixedPart
         }
         else{
         	d.coef<- NA
         	}
         	
        # add prediction error to the condition mean from the actual data
        g.conditional <- ghat + (g.unconditional.grid -  ghat.synthetic)
        # NOTE: sampling variablity for fixed part is built in too
        # because d.coef are estimated and included in the prediction. 
return(
       list( g.conditional = g.conditional, d.coef = d.coef) )
}
