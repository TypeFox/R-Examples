# Setting package options that need to be available on load.
.onAttach <- function (lib, pkg) 
{
  # Based on Hadley's template in devtools.
  # https://github.com/hadley/devtools/blob/master/R/zzz.r
  default.options <- options()
  managelocalrepo.options <- list(
    managelocalrepo.base = NA,
    managelocalrepo.submissions = NA
  )
  
  # Avoids overwriting options set by the user (perhaps using .First())
  to_set <- !(names(managelocalrepo.options) %in% names(default.options))
  if (any(to_set)) {
    options(managelocalrepo.options[to_set])
    packageStartupMessage('managelocalrepo: The following global options have ',
                          'been set:\n', 
                          paste(names(managelocalrepo.options), collapse='\n'))
  } 
}