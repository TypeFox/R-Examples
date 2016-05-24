
#' Create a new options manager.
#'
#' Set up a set of options with default values and retrieve a manager for it.
#'
#'
#' @section Details:
#' 
#' The function \code{options_manager} creates an option management function. The returned
#' function can be uset to set, get, or reset options. The only restriction of the package is
#' that the following words cannot be used as names for options:
#' 
#' \code{.__reset} \code{.__defaults}
#'
#' For more details and extensive examples see the vignette by copy-pasting this command:
#' 
#'    \code{vignette("settings", package = "settings")}
#'
#' @section Checking options:
#' Option values can be checked automatically by supplying the options manager
#' with a named list of functions (\code{.allowed}) that take an option value
#' and throw an error if it is out-of-range. The functions \code{\link{inlist}}
#' and \code{\link{inrange}} are convenience functions that create such checking
#' functions for you.
#'
#'
#' @param ... Comma separated \code{[name]=[value]} pairs. These will be the names and default values for your options manager.
#' @param .allowed list of named functions that check an option (see 'checking options') 
#'
#' @return A \code{function} that can be used as a custom options manager. It takes as arguments
#' a comma separated list of option names (\code{character}) to retrieve options or 
#' \code{[name]=[value]} pairs to set options.
#'   
#' @examples
#' # create an options register
#' my_options <- options_manager(foo=1,bar=2,baz='bob')
#' 
#' ### Retrieving options
#' my_options() # retrieve the full options list.
#' my_options('baz')
#' my_options('foo')
#' 
#' # When multiple options are retrieved, the result is a list
#' my_options('foo','baz')
#' 
#' ### Setting global options
#' my_options(foo=3,baz='pete')
#' my_options()
#' ### Reset options to default.
#' reset(my_options)
#' my_options()
#' 
#' ### Limit the possible values for an option.
#' my_options <- options_manager( fu="bar",.allowed = list(fu=inlist("foo","bar")) )
#'
#' @seealso 
#' 
#' Reset to default values: \code{\link{reset}}.
#' 
#' Retrieve default values: \code{\link{defaults}}
#' 
#' Create a local, possibly altered copy: \code{\link{clone_and_merge}}
#' 
#' @export
options_manager <- function(..., .allowed){
  stop_if_reserved(...)
  .defaults <- list(...)
  .op <- .defaults

  .al <- list()
  # default option checkers
  for ( v in names(.defaults)) .al[[v]] <- nolimit
  # set explicitly mentioned option checkers
  if (!missing(.allowed)) .al[names(.allowed)] <- .allowed
  # if you impose a limit on an option thats not defined: stop.
  if (!all(names(.al) %in% names(.op))  ){
    nm <- names(.al)[!names(.al) %in% names(.op)]
    stop(sprintf("Trying to set limits for undefined options %s\n",paste(nm,collapse=", ")))
  }
  # Check the default values against allowed ranges.
  vars <- names(.op)
  for (v in vars) .al[[v]](.defaults[[v]])
  
  # create the options manager function
  function(..., .__defaults=FALSE, .__reset=FALSE){
    L <- list(...)
    if (.__defaults) return(.defaults)
    if (.__reset){
      .op <<- .defaults
      return(invisible(.op))
    }
    # get all options
    if (length(L) == 0) return(.op)
    # set options:
    vars <- names(L)
    if ( !is.null(vars) && !any(vars == "") ){
      if (!all(vars %in% names(.defaults))){
        ii <- vars %in% names(.defaults)
        warning(sprintf("Ignoring options not defined in manager: %s"
            , paste(vars[!ii], collapse=", ") 
        ))
        vars <- vars[ii]
        L <- L[ii]
      }
      # check if values are allowed (only for occurring options).
      ii <- vars %in% names(.defaults)
      for ( v in vars[ii] ) .al[[v]](L[[v]])
      .op[vars] <<- L
      return(invisible(.op))
    }
    # get options
    if (is.null(vars)){
      vars <- unlist(L)
      return( if (length(vars)==1) .op[[vars]] else .op[vars] )
    }      
    stop("Illegal arguments")
  }
}


#' Option checkers
#'
#' These functions return a function that is used by the options manager internally
#' to check whether an option set by the user is allowed.
#'
#' @param ... comma-separated list of allowed values.
#' @param min minimum value (for numeric options)
#' @param max maximum value (for numeric options)
#' @seealso \code{\link{options_manager}} for examples.
#' @export
inlist <- function(...){
  .list <- unlist(list(...))
  function(x){
    if (!x %in% .list){
      stop(sprintf("Option value out of range. Allowed values are %s",paste(.list,collapse=",")),call.=FALSE)
    }
  }
}

#' @rdname inlist
#' @export
inrange <- function(min=-Inf,max=Inf){
  .range <- c(min=min, max=max)
  function(x){
    if( !is.numeric(x) || ( x > .range['max'] | x < .range['min']) ){
      stop(sprintf("Option value out of range. Allowed values are in [%g, %g]",.range['min'], .range['max'])
           ,call.=FALSE)
    }
  }
}

nolimit <- function(...) invisible(NULL) 




#' Create a local, altered copy of an options manager
#'
#' Local options management.
#'
#' @section Details:
#' This function creates a copy of the options manager \code{options}, with the same defaults.
#' However, the current settings may be altered by passing extra arguments. Its intended use
#' is to allow for easy merging of local options with global settings in a function call.
#' 
#' Some more examples can be found in the vignette: \code{vignette('settings',package='options')}.
#'
#' @param options A function as returned by \code{\link{options_manager}} or \code{clone_and_merge}.
#' @param ... Options to be merged, in the form of \code{[name]=[value]} pairs. 
#'
#' @return A option manager like \code{options}, with possibly different settings.
#'
#' @seealso \code{\link{options_manager}}, \code{\link{reset}}, \code{\link{defaults}}
#' 
#' @examples 
#' # Create global option manager.
#' opt <- options_manager(foo=1,bar='a')
#' 
#' # create an altered copy
#' loc_opt <- clone_and_merge(opt, foo=2)
#' 
#' # this has no effect on the 'global' version
#' opt()
#' # but the local version is different
#' loc_opt()
#' 
#' # we alter the global version and reset the local version
#' opt(foo=3)
#' reset(loc_opt)
#' opt()
#' loc_opt()
#' 
#' # create an options manager with some option values limited
#' opt <- options_manager(prob=0.5,y='foo',z=1,
#'   .allowed=list(
#'      prob = inrange(min=0,max=1)
#'      , y    = inlist("foo","bar")
#'    )
#'  )
#' # change an option
#' opt(prob=0.8)
#' opt("prob")
#' \dontrun{
#' # this gives an error
#' opt(prob=2)
#' }
#' 
#' 
#' @export 
clone_and_merge <- function(options,...){
  df <- options(.__defaults=TRUE)
  op <- options()
  f <- do.call(options_manager,df)
  do.call(f,op)
  f(...)
  f
}

#' Reset options to default values
#'
#' @param options An option manager, as returned by \code{\link{options_manager}} or \code{\link{clone_and_merge}}
#' 
#' @return The list of reset options, invisibly.
#' 
#' @seealso \code{\link{defaults}}
#' @export 
reset <- function(options) options(.__reset=TRUE)

#' Request default option values
#' 
#' @param options An option manager, as returned by \code{\link{options_manager}} or \code{\link{clone_and_merge}}
#' 
#' @return A \code{list}.
#' 
#' @seealso \code{\link{reset}}
#' @export
defaults <- function(options) options(.__defaults=TRUE)


#' Check for reserved option names.
#' 
#' Utility function checking for reserved names.
#' 
#' @section Details:
#' This is a utility function that checks if the keys of the key-value pairs
#' \code{...} contain reserved words. The reserved words are
#' 
#' \code{.__defaults}, \code{.__reserved}.
#' 
#' If reserved words are encountered in the input an error thrown.
#' The package vignette has examples of its use: 
#' 
#'    \code{vignette('settings',package='options')}
#' 
#' @param ... Comma-separated \code{[key]=[value]} pairs
#' 
#' @return \code{logical}, indicating if any of the keys was reserved (invisibly).
#' 
#' @seealso \code{\link{is_setting}}
#' 
#' @export
stop_if_reserved <- function(...){
  res <- c(".__defaults",".__reserved")
  out <- names(list(...)) %in% res
  if (any(out)){
    v <- paste(names(list(...))[out],collapse=", ")
    stop("Reserved word used as option name: ",v)
  }
  invisible(out)
}

#' Find out if we're setting or getting
#' 
#' Utility function checking if we're setting or getting.
#' 
#' @param ... \code{[key]=[value]} pairs of options
#' @return \code{logical}, \code{TRUE} if \code{...} represents set-options, \code{FALSE} if
#'  \code{...} represents get-options. An error is thrown if it cannot be determined.
#' 
#' 
#' @seealso \code{\link{stop_if_reserved}}
#' 
#' @export 
is_setting <- function(...){
  L <- list(...)
  nm <- names(L)
  set <- !is.null(nm) && !any(nm=="")
  get <- length(L)  == 0 | (length(L) > 0 & is.null(nm)) 
  stopifnot(set | get)
  set
}


