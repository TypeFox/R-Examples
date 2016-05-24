#' Load a list of libraries, reporting to the Rattle Log
#'
#' Only load the package if not already loaded. If already loaded then
#' we don't return the name from the function.
#'
#' @param l Vector of pairs, "package name" "used function".
#' @return returns list of packages that get loaded.
#' @rdname loadLibs
loadLibs <- function(l)
{
  odd   <- seq(1, length(l), 2)
  lname <- l[odd]
  even  <- seq(2, length(l), 2)
  lfun  <- l[even]
  libs  <- NULL
  for (i in 1:length(odd))
  {
    appendLog(packageProvides(lname[i], lfun[i]), sprintf("library(%s)", lname[i]))
    if (!sprintf("package:%s", lname[i]) %in% search())
    {
      suppressPackageStartupMessages(library(lname[i], character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE))
      libs <- c(libs, lname[i])
    }
  }
  return(libs)
}

