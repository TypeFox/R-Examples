##' Run the unit tests
##'
##' Run the unit tests attached to the functions via \link[svUnit]{svUnit} 
##' @return invisibly \code{TRUE} if the tests pass, \code{NA} if \link[svUnit]{svUnit} is not
##' available. Stops if errors are encountered.
##' @author Claudia Beleites
##' @seealso  \link[svUnit]{svUnit} 
##' @export 
arrayhelpers.unittest <- function (){
  if (! require ("svUnit", quietly = TRUE)){
    warning ("svUnit required to run the unit tests.")
    return (NA)
  }

  tests <- unlist (eapply (env = getNamespace ("arrayhelpers"), FUN = is.test, all.names = TRUE))
  tests <- names (tests [tests])
  tests <- sapply (tests, get, envir = getNamespace ("arrayhelpers"))

  clearLog ()
  warnlevel <- options()$warn
  options (warn = 0)
  for (t in seq_along (tests))
    runTest (tests [[t]], names (tests) [t])

  options (warn = warnlevel)
  if (interactive ())
    print (stats (Log()))
  else
    print (stats (Log ())[,c ("kind", "msg")])

  errorLog (summarize = FALSE)
  invisible (TRUE)
}


