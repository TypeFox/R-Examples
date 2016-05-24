#' @title Various environment manipulation functions.
#' 
#' @description The computation of \emph{adaptive} power spectral density estimates 
#' requires bookkeeping and non-destructive manipulation of variables.  
#' The functions here are mainly convenience wrappers
#' designed to maintain variable separation from the 
#' \code{.GlobalEnv} environment so that no innocent variable is destroyed in
#' the process of iteratively computing spectra.
#' \strong{The user should generally not be using the \emph{setters} even though
#' all functions exist in the namespace.}
#' 
#' \code{get_psd_env_pointer} is a convenience wrapper to get the environment pointer.
#' 
#' \code{get_psd_env_name} is a convenience wrapper to get the environment name.
#'
#' \code{psd_envRefresh} will clear any variables in the enviroment and reset the initialization stamp.
#' 
#' \code{psd_envClear} clears the contents of the environment.
#' 
#' \code{psd_envStatus} returns a list of some information regarding
#' the status of the environment.
#' 
#' \code{psd_envList} returns a listing of any assignments.
#' 
#' \code{psd_envGet} returns the value of \code{variable}.
#' 
#' \code{psd_envAssign} assigns \code{value} to \code{variable}, but does not return it.
#' 
#' \code{psd_envAssignGet} both assigns and returns a value.
#' 
#' \code{update_adapt_history} updates the adaptive estimation history list.
#' 
#' \code{new_adapt_history} initializes a nested-list object to store the 
#' data from each iteration.
#' 
#' @rdname psd-environment
#' @name psd-environment
#' @docType methods
#' @section Defaults and Initialization:
#' One can use \code{get_psd_env_pointer()} and \code{get_psd_env_name()} to access the 
#' pointer and name of the environment, if
#' needed.
#'
#' \code{psd_envRefresh} should be used when
#' a fresh environment is desired: typically only if, for example, \code{\link{psdcore}} is 
#' used rather than \code{\link{pspectrum}}.
#'
#' @section Assigning and Retieving:
#' \code{psd_envAssign} and \code{psd_envGet} perform the assignments and retrieval
#' of objects in the environment.  A convenience function, \code{psd_envAssignGet},
#' is included so that both assignment and retrieval may be performed at the same
#' time.  This ensures the assignment has succeeded, and the returned value is
#' not from some other frame.
#'
#' @section Getters and Setters:
#' The functions here can be classified whether the get, or set variables in the
#' environment; some do both.  
#' Others make no modifications to the environment.
#' 
#' \subsection{Getter}{
#' \itemize{
#' \item{\code{get_adapt_history}}{}
#' \item{\code{get_psd_env_name}}{}
#' \item{\code{get_psd_env_pointer}}{}
#' \item{\code{psd_envGet}}{}
#' \item{\code{psd_envList}}{}
#' \item{\code{psd_envStatus}}{}
#' }
#' }
#'
#' \subsection{Setter}{
#' \itemize{
#' \item{\code{new_adapt_history}}{}
#' \item{\code{psd_envAssign}}{}
#' }
#' }
#' 
#' \subsection{Getter and Setter}{
#' \itemize{
#' \item{\code{psd_envAssignGet}}{}
#' \item{\code{psd_envClear}}{}
#' \item{\code{psd_envRefresh}}{}
#' \item{\code{update_adapt_history}}{}
#' }
#' }
#'
#' @param verbose logical; should messages be given?
#' @param variable character; the name of the variable to get or assign
#' @param value character; the name of the variable to assign
#' @param adapt_stages scalar; The number of adaptive iterations to save (excluding pilot spectrum).
#' @param stage scalar; the current stage of the adaptive estimation procedure
#' @param ntap vector; the tapers
#' @param PSD vector or object with class \code{'spec'}; the power spectral density estimates
#' @param freq vector; the frequencies
#' @param ... additional arguments
#' 
#' @seealso \code{\link{psd-utilities}}, \code{\link{pspectrum}}
#' @example inst/Examples/rdex_psdenv.R
NULL

#' @rdname psd-environment
#' @export
get_psd_env_pointer <- function() .psdEnv

#' @rdname psd-environment
#' @export
get_psd_env_name <- function() .psdEnvName

#' @rdname psd-environment
#' @export
psd_envRefresh <- function(verbose=TRUE) {
  # env params
  envname <- get_psd_env_name()
  # rm all in envir
  psd_envClear()
  psd_envAssign("init", sprintf("refreshed at %s", Sys.time()))
  if (verbose) {
    message(sprintf("\tenvironment  ** %s **  refreshed", envname))
  }
  return(invisible(psd_envStatus()))
}

#' @note \code{psd_envClear} does \emph{not} remove the environment--simply the assignments within it.
#' @rdname psd-environment
#' @export
psd_envClear <- function(){
  ENV <- get_psd_env_pointer()
  listing <- psd_envList()
  rm(list=listing, envir=ENV)
}

#' @rdname psd-environment
#' @export
psd_envStatus <- function(){
  envname <- get_psd_env_name()
  envir <- get_psd_env_pointer()
  return(list(env_name=envname, 
              env_pointer=envir, 
              env_is_env=is.environment(envir), 
              listing=psd_envList(),
              env_init=psd_envGet("init"),
              env_status_stamp=Sys.time() ))
  
}

#' @rdname psd-environment
#' @export
psd_envList <- function(){
  ## return listing of envir::variable
  ENV <- get_psd_env_pointer()
  ls(envir=ENV, all.names=TRUE)
}

#' @rdname psd-environment
#' @export
psd_envGet <- function(variable){
  ## return contents on envir::variable
  ENV <- get_psd_env_pointer()
  if (!exists(variable, envir=ENV)){
    warning(sprintf("Variable  '%s'  not found!  Inspect results of  psd_envList()", variable))
    return(NULL)
  } else {
    return(get(variable, envir=ENV))
  }
}

#' @rdname psd-environment
#' @export
psd_envAssign <- function(variable, value){
  ## set contents of envir::variable to value
  ENV <- get_psd_env_pointer()
  assign(variable, value, envir=ENV)
}

#' @rdname psd-environment
#' @export
psd_envAssignGet <- function(variable, value){
  ## set contents of envir::variable to value
  psd_envAssign(variable, value)
  psd_envGet(variable)
}

#' @export
#' @rdname psd-environment
get_adapt_history <- function() psd_envGet("histlist")

#' @export
#' @rdname psd-environment
last_psd <- function(){
  fp <- psd_envGet("final_psd")
  if (is.null(fp)) stop('last psd not available')
  return(fp)
}

#' @rdname psd-environment
#' @export
update_adapt_history <- function(PSD, stage, ...) UseMethod("update_adapt_history")

#' @rdname psd-environment
#' @export
update_adapt_history.spec <- function(PSD, stage, ...){
  p. <- PSD[['spec']]
  s. <- stage
  t. <- as.tapers(PSD[['taper']])
  f. <- PSD[['freq']]
  update_adapt_history(p., s., t., f., ...)
}

#' @rdname psd-environment
#' @export
update_adapt_history.default <- function(PSD, stage, ntap=NA, freq=NULL, ...){
  stopifnot(stage >= 0)
  histlist <- get_adapt_history()
  # stage == 0 <--> index == 1
  stg_ind <- stage + 1
  if (!is.null(freq)) histlist[['freq']] <- freq
  histlist[['stg_kopt']][[stg_ind]] <- ntap
  histlist[['stg_psd']][[stg_ind]] <- PSD
  if (is.null(histlist[['freq']]) & stage > 0) warning("frequencies are absent despite non-pilot stage update")
  psd_envAssignGet("histlist", histlist)
}


#' @rdname psd-environment
#'
#' @section Adaptive History:
#' The list object for historical adapt-data may be accessed with \code{\link{get_adapt_history}}.
#' The top names of the returned list are
#' \describe{
#' \item{\code{stg_kopt}}{Sequential taper vectors.}
#' \item{\code{stg_psd}}{Sequential power spectral density vectors.}
#' \item{\code{freq}}{The frequencies for each set of \code{stg_kopt} and \code{stg_psd}.}
#' }
#' @export
new_adapt_history <- function(adapt_stages){
  stopifnot(length(adapt_stages)==1)
  histlist <- vector("list", 3) # freq, list-tap, list-psd
  names(histlist) <- c("freq", "stg_kopt", "stg_psd")
  num_pos <- 1 + adapt_stages # pilot + adapts
  histlist[[2]] <- histlist[[3]] <- vector("list", num_pos)
  psd_envAssignGet("histlist", histlist)
}
