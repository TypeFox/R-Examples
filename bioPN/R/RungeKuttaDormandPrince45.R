## -*- ess-indent-level: 2; ess-basic-offset: 2; tab-width: 8 -*-
##
## Copyright (C) 2009-2014 Roberto Bertolusso and Marek Kimmel
##
## This file is part of bioPN.
##
## bioPN is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## bioPN is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with bioPN. If not, see <http://www.gnu.org/licenses/>.


RungeKuttaDormandPrince45 <- function(model, timep, delta=1, ect=1e-9) {
  model$slow=rep(0,dim(model$pre)[1])

  return(HaseltineRawlings(model, timep, delta, runs=1, ect=ect))
}
