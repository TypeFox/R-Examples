# DESP/R/scsSOCP.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (version 3) as published by
#  the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

scsSOCP <- 
function(model){
  # solve an SOCP using SCS (direct method) :

  # settings
  lst = list("alpha"=-1, "rho_x"=-1, "max_iters"=-1, "scale"=-1, "eps"=-1, "cg_rate"=-1, "verbose"=-1, "normalize"=-1, "warm_start"=-1)
  if (is.null(model$settings)) model$settings <- lst
  else {
    nSett <- names(lst)
    model$settings <- mapply(function(a,b) ifelse(is.null(b),a,b),lst, model$settings[nSett], SIMPLIFY=FALSE)
  }

  # call SCS
  .Call("scs_SOCP_solve",as.double(model$Ax),as.integer(model$Ai),as.integer(model$Ap),as.integer(model$Am),as.integer(model$An),as.double(model$b),as.double(model$c),as.integer(model$Kf),as.integer(model$Kl),as.integer(model$Kq),as.integer(model$Kqsize),as.double(model$settings$alpha),as.double(model$settings$rho_x),as.integer(model$settings$max_iters),as.double(model$settings$scale),as.double(model$settings$eps),as.double(model$settings$cg_rate),as.integer(model$settings$verbose),as.integer(model$settings$normalize),as.integer(model$settings$warm_start))
}
