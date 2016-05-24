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


GillespieOptimDirect <- function(model, timep, delta=1, runs=1) {
  pre <- model$pre
  post <- model$post
  h <- model$h
  M <- model$M

  storage.mode(pre) <- storage.mode(post) <- storage.mode(runs) <- "integer"
  storage.mode(delta) <- storage.mode(M) <- "double"

  if (is.null(place <- model$place)) {
    place <- paste("P",1:ncol(pre),sep="")
  }
  if (is.null(transition <- model$transition)) {
    transition <- paste("T",1:nrow(pre),sep="")
  }

  .Call("GillespieOptimDirect", pre, post, h, M, timep, delta, runs, place, transition, parent.frame())
}
