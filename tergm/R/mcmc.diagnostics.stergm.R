#  File R/mcmc.diagnostics.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
mcmc.diagnostics.stergm <- function(object, 
                                    center=TRUE,
                                    curved=TRUE,
                                    vars.per.page=3, ...){
  if(!is.null(object$formation.fit$sample)){
    cat("\n==========================\n")
    cat("Formation fit diagnostics\n")
    cat("==========================\n\n")
    mcmc.diagnostics.ergm(object$formation.fit, center=center, curved=curved, vars.per.page=vars.per.page, ...)
  }
  if(!is.null(object$dissolution.fit$sample)){
    cat("\n==========================\n")
    cat("Dissolution fit diagnostics\n")
    cat("==========================\n\n")
    mcmc.diagnostics.ergm(object$dissolution.fit, center=center, curved=curved, vars.per.page=vars.per.page, ...)
  }
  if(!is.null(object$sample)){
    cat("\n==========================\n")
    cat("EGMME diagnostics\n")
    cat("==========================\n\n")
    mcmc.diagnostics.ergm(object, center=center, curved=curved, vars.per.page=vars.per.page, ...)
  }
}
