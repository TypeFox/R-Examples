## -*- ess-indent-level: 2; ess-basic-offset: 2; tab-width: 8 -*-
##
## Copyright (C) 2009-2014 Roberto Bertolusso and Marek Kimmel
##
## This file is part of sbioPN.
##
## sbioPN is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## sbioPN is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with sbioPN. If not, see <http://www.gnu.org/licenses/>.

sHaseltineRawlings <- function(model, timep, delta=1, runs=1, ect=1e-9) {
  storage.mode(runs) <- storage.mode(burnRnd) <- "integer"
  storage.mode(delta) <- storage.mode(ect) <- "double"

  if (is.null(place <- model$place)) {
    place <- paste("P",1:13,sep="")
  }
  if (is.null(transition <- model$transition)) {
    transition <- paste("T",1:13,sep="")
  }

  model <- list(VoxelFamily=list(model))
  .Call("sHaseltineRawlings", model, timep, delta, runs, place, transition, ect, 0, parent.frame())
}
