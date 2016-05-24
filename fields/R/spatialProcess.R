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
spatialProcess <- function(x, y, cov.function = "stationary.cov", 
	cov.args = list(Covariance = "Matern", smoothness = 1),
	 ngrid=10, theta.grid = NULL, ...) {
	MLEfit <- MLESpatialProcess(x, y, cov.function = cov.function, 
		cov.args = cov.args, ngrid=ngrid, theta.grid = theta.grid,
		 ...)
# now fit spatial model with MLE for theta (range parameter)
# reestimate the other parameters for simplicity to get the complete Krig
# object.		
	obj <- Krig(x, y, cov.function = cov.function, cov.args = cov.args, 
		theta = MLEfit$pars[1],  
		method = "REML", give.warnings=TRUE, 
		...)
	obj <- c(obj, MLEfit)
	obj$theta.MLE<- MLEfit$pars[1]
# replace call with this top level one
    obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","Krig")
 
	return(obj)
}
