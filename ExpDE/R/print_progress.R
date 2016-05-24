#' Print progress of DE
#'
#' Echoes the progress of DE to the terminal
#'
#' @section Parameters:
#' This routine accesses all variables defined in the calling environment using
#' \code{parent.frame()}, so it does not require any explicit input parameters.
#' However, the calling environment must contain:
#' \itemize{
#'    \item \code{showpars}: list containing parameters that control the printed
#'                           output of \code{moead()}. Parameter \code{showpars}
#'                           can have the following fields:
#'    \itemize{
#'      \item \code{$show.iters = c("dots", "numbers", "none")}: type of output.
#'                                Defaults to \code{"numbers"}.
#'      \item \code{$showevery}: positive integer that determines how frequently
#'                               the routine echoes something to the terminal.
#'                               Defaults to \code{1}.
#'    }
#'    \item \code{iters()} : counter function that registers the iteration
#'                           number
#' }
#'
#' @export
#' 
print_progress <- function(){
  
  # Get access to the variables in the calling environment
  env  <- parent.frame()
  
  # ========== Error catching and default value definitions
  #tmp <- assert_that(has_name(env, "showpars"))
  pars <- env$showpars
  
  if(!any("show.iters" == names(pars))) pars$show.iters <- "numbers"
  if(!any("showevery" == names(pars)))  pars$showevery  <- 1
  if(!any("show.plot" == names(pars)))  pars$show.plot  <- FALSE
  
  #   tmp <- assert_that(
  #     is.count(pars$showevery),
  #     any(pars$show.iters == c("dots", "numbers", "none")),
  #     length(pars$show.iters) == 1,
  #     length(pars$showevery) == 1,
  #     is.flag(pars$show.plot))
  # ==========
  
  if (pars$show.iters != "none"){
    if (env$t == 1) cat("\nExpDE running: ")
    if (env$t %% pars$showevery == 0){
      if (pars$show.iters == "dots") cat(".")
      if (pars$show.iters == "numbers") cat("\nIteration: ", env$t)
    }
  }
}