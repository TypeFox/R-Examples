#' DYM
#' 
#' You might mistype an object name.
#' The package suggests the correct spell of the object you meant.
#' 
#' @name DYM
#' @docType package
NULL

#' If the function is called after an error of 'object not found',
#' the function tries to tell you the name of the correct name that you meant.
#' 
#' @param threshold
#'    The maximum distance between the misspell (\code{x}) and the correct answer (in \code{name}).
#' @param max_n
#'    An integer limiting the number of results.  Passed to \code{\link[utils]{head}}.
#' @param ignoreCase
#'    A logical value indicating whether differences in case should be ignored when matching.  Passed to \code{\link[utils]{adist}}.
#' @examples 
#' \dontrun{
#' options(error = DYM::DYM())
#' logg(10)
#' 
#' # For fewer or more suggestions, change threshold, max_n and ignoreCase
#' options(error = DYM::DYM(threshold = 3, max_n = 25, ignoreCase = TRUE))
#' logg(10)
#' }
#' @export
DYM <- function(threshold = 2, max_n = 10, ignoreCase = FALSE)
{
   threshold <- force(threshold)
   max_n <- force(max_n)
   ignoreCase <- force(ignoreCase)
   function()
   {
      missingVariable <- getMissingVariable()
      if (!is.na(missingVariable)) {
         envir <- attr(missingVariable, "package")
         if (is.null(envir)) {
            envir <- .GlobalEnv
         }
         availableVariables <- getNames(class(missingVariable), envir)
         names <- findSimilarName(missingVariable, availableVariables, threshold, max_n, ignoreCase)
         if (length(names) > 0L) {
            message <- ngettext(length(names), "Did you mean: %s", "Did you mean: [%s]", domain="R-DYM")
            hints <- sapply(names, sprintf, fmt="'%s'")  # sQuote might be better
            message(sprintf(message, paste(hints, collapse=", ")))
         }
      }
   }
}
