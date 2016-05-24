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

ssa.check.method <- function(x0,a,nu,method,tau,f) {

  # Check the consistency of the system dimensions, i.e. number of rows and 
  # columns in the state-change matrix and the number of elements in the initial 
  # state vector and the vector of propensity functions  
  if ((length(a)/dim(nu)[2]) != (length(x0)/dim(nu)[1])) 
    stop("inconsistent system dimensions (unequal 'nu' tessallation)")
  if (((length(a)%%dim(nu)[2])>0) || ((length(x0)%%dim(nu)[1])>0)) 
  stop("inconsistent system dimensions (fractional tessallation)")

  # For the ETL method tau>0
  if ((method=="ETL") & (!(tau>0))) stop("ETL method requires tau>0") 

  # Check that f (used in the BTL method) is >1 
  if (method=="BTL" & f<=1) stop("f has to be >1") 
}