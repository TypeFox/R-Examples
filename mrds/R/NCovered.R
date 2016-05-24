#' Compute estimated abundance in covered (sampled) region
#'
#' Generic function that computes abundance within the covered region.  It
#' calls method (class) specific functions for the computation.
#'
#' @aliases NCovered NCovered.ds NCovered.io NCovered.io.fi NCovered.trial
#'   NCovered.trial.fi NCovered.rem NCovered.rem.fi
#' @param par parameter values (used when computing derivatives wrt parameter
#'   uncertainty); if NULL parameter values in \code{model} are used
#' @param model ddf model object
#' @param group if TRUE computes group abundance and if FALSE individual
#'   abundance
#' @return abundance estimate
#' @author Jeff Laake
#' @keywords utility
NCovered <- function(par,model=NULL,group=TRUE){
  # Functions Used: NCovered.ds, NCovered.io, NCovered.io.fi, NCovered.trial,
  #                 NCovered.trial.fi, NCovered.rem, NCovered.rem.fi

  if(is.null(model)){
    model <- par
    par <- NULL
  }

  # call method specific function
  result <- switch(model$method,
         ds=NCovered.ds(par=par,model=model,group=group),
         io=NCovered.io(par=par,model=model,group=group),
         io.fi=NCovered.io.fi(par=par,model=model,group=group),
         trial=NCovered.trial(par=par,model=model,group=group),
         trial.fi=NCovered.trial.fi(par=par,model=model,group=group),
         rem=NCovered.rem(par=par,model=model,group=group),
         rem.fi=NCovered.rem.fi(par=par,model=model,group=group))
  return(result)
}
