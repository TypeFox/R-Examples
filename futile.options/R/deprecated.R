# This is here for backwards compatibility
options.manager <- function(option.name, defaults=NULL)
{
  OptionsManager(option.name, defaults)
}

# This is here for backwards compatibility
reset.options <- function(option.name, ...) resetOptions(option.name, ...)

