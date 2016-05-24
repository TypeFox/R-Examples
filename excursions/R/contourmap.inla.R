## excursions.inla.R
##
##   Copyright (C) 2012, 2013, 2014 David Bolin, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.


contourmap.inla <- function(result.inla,
                            stack,
                            name=NULL,
                            tag=NULL,
                            method="EB",
                            ind,...)
{
  if (!requireNamespace("INLA", quietly=TRUE))
    stop('This function requires the INLA package (see www.r-inla.org/download)')
  if(missing(result.inla))
    stop('Must supply INLA result object')

  if(method != "EB")
    stop("Currently only EB method is implemented")

  ind.stack <- inla.output.indices(result.inla, name=name, stack=stack, tag=tag)
  n.out <- length(ind.stack)
  ind.int <- seq_len(n.out)
  if(!missing(ind) && !is.null(ind)){
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
  }
  ind = ind.stack

  for(i in 1:result.inla$misc$configs$nconfig){
    config = private.get.config(result.inla,i)
    if(config$lp == 0)
      break
  }
  cm <- contourmap(mu=config$mu,Q = config$Q, ind=ind,...)

  set.out = rep(NA,n.out)
  if(!is.null(cm$E)){
    set.out[ind.int] = cm$E[ind];
    cm$E = set.out
  }

  if(!is.null(cm$G)) {
    set.out[ind.int] = cm$G[ind];
    cm$G = set.out
  }

  if(!is.null(cm$F)){
    F0 = cm$F[ind]; F0[is.na(F0)] = 0
    cm$P0 = mean(F0)
    set.out[ind.int] = cm$F[ind];
    cm$F = set.out
  }
  if(!is.null(cm$M))
    set.out[ind.int] = cm$M[ind]; cm$M = set.out

  if(!is.null(cm$rho))
    set.out[ind.int] = cm$rho[ind]; cm$rho = set.out

  if(!is.null(cm$map))
    set.out[ind.int] = cm$map[ind]; cm$map = set.out

  cm$meta$ind <- ind.int
  cm$meta$call = match.call()
  return(cm)
}
