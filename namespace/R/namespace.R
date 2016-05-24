# Returns the namespace registry
#' @useDynLib devtools nsreg
.getNameSpaceRegistry <- function() {
  .Call("nsreg", "namespace")
}

# Register a namespace
registerNamespace <- function(name = NULL, env = NULL) {
  # Be careful about what we allow
  if (!is.character(name) || name == "" || length(name) != 1)
    stop("'name' must be a non-empty character string.")

  if (!is.environment(env))
    stop("'env' must be an environment.")

  if (name %in% loadedNamespaces())
    stop("Namespace ", name, " is already registered.")

  # Add the environment to the registry
  nsr <- .getNameSpaceRegistry()
  nsr[[name]] <- env

  env
}


# Unregister a namespace - should be used only if unloadNamespace()
# fails for some reason
unregisterNamespace <- function(name = NULL) {
  # Be careful about what we allow
  if (!is.character(name) || name == "" || length(name) != 1)
    stop("'name' must be a non-empty character string.")

  if (!(name %in% loadedNamespaces()))
    stop(name, " is not a registered namespace.")

  # Remove the item from the registry
  rm(name, envir=.getNameSpaceRegistry())
  invisible()
}

# This is similar to getNamespace(), except that getNamespace will load
# the namespace if it's not already loaded. This function will not.
# In R 2.16, a function called .getNamespace() will have the same effect
# and this will no longer be necessary.
if(exists(".getNamespace", where="package:base", mode="function")) {
  getRegisteredNamespace <- function(name) base::.getNamespace(name)
} else {
  getRegisteredNamespace <- function(name)
  {
    ## Sometimes we'll be passed something like as.name(name), so make sure
    ## it's a string for comparison
    if (!(as.character(name) %in% loadedNamespaces()))
      return(NULL)
    else
      return(getNamespace(name))
  }
}
