#' install.package and require
#' 
#' install and load a package. If a package is not available, it is installed before being loaded
#' 
#' @aliases require2 library2
#' @return Cats help instruction.
#' @note Passing a vector with packages will work, but give some warnings.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2014
#' @seealso \code{\link{install.packages}}, \code{\link{require}}
#' @keywords package
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Excluded fom CRAN checks. Package installation on server is unnecessary.
#' require2(ada)
#' }
#' 
#' @param name Name of the package(s). Can be qouted, must not.
#' @param \dots Arguments passed to \code{\link{install.packages}} like \code{lib}, \code{repos} etc.
#' 
require2 <- function(
name,
...)
{
name <- as.character(substitute(name))
for(i in 1:length(name))
{
if(!requireNamespace(name[i], quietly=TRUE))
   {
   install.packages(name, ...)
   library(name[i], character.only=TRUE)
   }
}
for(i in 1:length(name))
  message(paste0('-------------------------\nhelp(package="', name[i],
            '")\n-------------------------\n'))
}

#' @export
library2 <- require2
