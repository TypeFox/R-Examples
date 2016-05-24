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

ssa.check.args <- function(x0,a,nu,tf,method,tau,f,epsilon,nc,hor,dtf,nd,ignoreNegativeState,consoleInterval,censusInterval,verbose) {

  # Do some basic check of the argument types
  if (!is.numeric(x0))              stop("'x0' is not numeric")
  if (!is.character(a))             stop("'a' is not of character type")
  if (!is.numeric(nu))              stop("'nu' is not numeric")
  if (!is.numeric(tf))              stop("'tf' is not numeric")
  if (!is.character(method))        stop("'method' is not of character type")
  if (!is.numeric(tau))             stop("'tau' is not numeric")
  if (!is.numeric(f))               stop("'f' is not numeric")
  if (!is.numeric(epsilon))         stop("'epsilon' is not numeric")
  if (!is.numeric(nc))              stop("'nc' is not numeric")
  if (!is.numeric(hor))             stop("'hor' is not numeric")
  if (!is.numeric(dtf))             stop("'dtf' is not numeric")
  if (!is.numeric(nd))              stop("'nd' is not numeric")
  if (!is.numeric(consoleInterval)) stop("'consoleInterval' is not numeric")
  if (!is.numeric(censusInterval))  stop("'censusInterval' is not numeric")
  if ((ignoreNegativeState != TRUE) & (ignoreNegativeState != FALSE)) 
    stop("'ignoreNegativeState' is not boolean")
  if ((verbose != TRUE) & (verbose != FALSE)) stop("'verbose' is not boolean")
  if (is.null(names(x0))) stop("'x0' is missing element names")
}
