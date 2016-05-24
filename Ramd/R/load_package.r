#' Load a package and install it if necessary
#' 
#' @name load_package
#' @param name Name of package
#' @examples
#' \dontrun{
#' load_package('glmnet')
#' }
load_package <- function(name) {
  if (!require(name, character.only = TRUE)) {
    install.packages(name, dependencies = TRUE)
    if (!require(name, character.only = TRUE))
      stop(paste('Package', name, 'not found'))
  }
}
