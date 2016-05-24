# Passing the update argument to the resulting options manager will update
# the provided keys and values. Note that update[[1]] are the keys and 
# update[[2]] are the values.
# Examples of using options.manager
#  log.options <- options.manager('log.options', defaults=list(logger='ROOT'))
#  log.options(a=123, b=6234)
#  log.options()
#  log.options(a=123, b=6234)
#  resetOptions(log.options, c=29)
#  log.options()
# Generates a function to retrieve options for a given name
OptionsManager <- function(option.name, defaults=list())
{
  # Define a function to update options
  up <- function(os)
  {
    my.options <- list()
    my.options[[option.name]] <- os
    options(my.options)
  }

  up(list())

  function(..., simplify=FALSE, update=list())
  {
    os <- getOption(option.name)
    if (length(os) < 1)
    {
      os <- defaults
      up(os)
    }

    args <- list(...)
    if (length(update) > 0)
      invisible(updateOptions(option.name, update[[1]], update[[2]]))
    if (length(args) == 0) return(os)

    # Getter
    if (any(is.null(names(args))))
    {
      ns <- sapply(args, '[')
      ns <- ns[ns %in% names(os)]
      if (length(ns) == 0) return(NULL)
      if (length(ns) == 1) return(os[[ns]])
      return(sapply(os[ns], '[', simplify=simplify))
    }

    # Setter
    for (x in names(args)) os[[x]] <- args[[x]]
    up(os)
    invisible()
  }
}

# Reset options for a given option set
resetOptions <- function(option.name, ...) UseMethod('resetOptions')
resetOptions.default <- function(option.name, ...)
  resetOptions.character(deparse(substitute(option.name)), ...)

resetOptions.character <- function(option.name, ...)
{
  my.options <- list()
  my.options[[option.name]] <- NA
  options(my.options)

  args <- list(...)
  if (length(args) > 0)
  {
    ks <- names(args)
    vs <- sapply(args, '[')
    kvs <- paste(ks,vs, sep='=')
    line <- paste(kvs, collapse=',')

    exp <- parse(text=paste('new.options <- ',option.name,'(',line,')',sep=''))
    eval(exp)
  }
  invisible()
}

# Update a specific option in an option set (a generated function)
# This is primarily used when dynamic creation of options variables are needed.
updateOptions <- function(option.name, ...) UseMethod('updateOptions')
updateOptions.default <- function(option.name, ...)
  updateOptions.character(deparse(substitute(option.name)), ...)

updateOptions.character <- function(option.name, key, value, ...)
{
  os <- getOption(option.name)
  if (is.null(os))
    stop(paste("Cannot update non-existent options:",option.name))

  if (length(key) == 1)
    os[[key]] <- value
  else
    for (idx in 1:length(key)) os[[key[idx]]] <- value[idx]
  my.options <- list()
  my.options[[option.name]] <- os
  options(my.options)

  invisible()
}

