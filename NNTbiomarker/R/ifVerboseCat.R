#' setVerboseCatOption
#'
#' Allows user to toggle on and off printing messages on a per-function basis.
#' Should be usable in other packages, but not by importing.
#'
#' @aliases ifVerboseCat
#' @param fname Name of the function to control.
#' @param value Boolean value: should this function print out messages?
#' @return The new value of the namespace option for fname ifVerboseCat
#'
setVerboseCatOption = function(fname, value) {
  namespaceName = getNamespaceName(environment(ifVerboseCat))
  opts = options()
  opts[[namespaceName]] [fname] = value
  options(opts)
  opts[[namespaceName]]
}


ifVerboseCat = function(...){
  # Because of this line, you can copy this code into another package project
  # without change, but you can't import this function.
  namespaceName = getNamespaceName(environment(ifVerboseCat))
  optionName = paste0(namespaceName, ".ifverboseCat")
  ### Initialize the option set.
  if(is.null(options(namespaceName)[[1]])) {
    setVerboseCatOption(fname, logical(0))
  }
  ### Get name of calling function.
  fname = try(as.character(parse(text=sys.call(-1)[1]))[1] )
  if(class(fname) == "try-error") return(invisible(NULL))
  ### If option not yet set, set it to TRUE.
  if(is.na(options(namespaceName)[[1]][fname])) {
    setVerboseCatOption(fname, TRUE)
  }
  ### If TRUE, print out the requested line of text.
  if(options(namespaceName)[[1]] [fname])
    catn(fname, ": ", ...)
  invisible(NULL)
}
