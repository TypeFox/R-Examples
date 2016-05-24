# Copyright 2007, 2008, 2010 Mario Pineda-Krch.
#
# This file is part of the R package GillespieSSA.
#
# GillespieSSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GillespieSSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GillespieSSA.  If not, see <http://www.gnu.org/licenses/>.

`ssa` <-
function(         x0 = stop("undefined 'x0'"), 
                   a = stop("undefined 'a'"),
                  nu = stop("undefined 'nu'"),
               parms = NULL, 
                  tf = stop("undefined 'tf'"), 
              method = "D",
             simName = "",
                 tau = 0.3,
                   f = 10,
             epsilon = 0.03,
                  nc = 10,
                 hor = NaN,
                 dtf = 10,
                  nd = 100,
 ignoreNegativeState = TRUE,
     consoleInterval = 0,
      censusInterval = 0,
             verbose = FALSE,
         maxWallTime = Inf
             ) 
{ # End of function arguments
 
  ssa.check.args(x0,a,nu,tf,method,tau,f,epsilon,nc,hor,dtf,nd,ignoreNegativeState,consoleInterval,censusInterval,verbose)
  
  # Convert lower case method names to upper case (undocumented featurette)
  if (method=="d")   method <- "D"
  if (method=="etl") method <- "ETL"
  if (method=="btl") method <- "BTL"
  if (method=="otl") method <- "OTL"

  ssa.check.method(x0,a,nu,method,tau,f)
  
  # Is the system nu-tiled along the diagonal?
  if ((length(a)/dim(nu)[2]>1) && (length(x0)/dim(nu)[1])>1){
    if (method=="D")   method <- "D.diag" 
    if (method=="ETL") method <- "ETL.diag"
    if (method=="BTL") method <- "BTL.diag"
    if (method=="OTL") method <- "OTL.diag"
  }
  
  # Take a snapshot of all the options so they can be saved later
  args <- list(    x0 = x0, 
                    a = a, 
                   nu = nu, 
                parms = parms, 
                   tf = tf, 
               method = method, 
                  tau = tau, 
                    f = f, 
              epsilon = epsilon, 
                   nc = nc, 
                  hor = hor, 
                  dtf = dtf, 
                   nd = nd, 
  ignoreNegativeState = ignoreNegativeState, 
      consoleInterval = consoleInterval, 
       censusInterval = censusInterval,
              verbose = verbose,
              simName = simName)
  
  # Run the simulation
  out.rxn <- ssa.run(x0,a,nu,parms,tf,method,tau,f,epsilon,nc,hor,dtf,nd,ignoreNegativeState,consoleInterval,censusInterval,verbose,maxWallTime) 

  # Wrap up the simulation
  out.summary <- ssa.terminate(args,out.rxn,tf,method,maxWallTime,verbose)
  return(out.summary)  
}
