##
## Replaces the (defunct) addLast() function.
##
lastAdd <- function( fun )
  {
    if (!is.function(fun)) stop("fun must be a function")
    if(!exists(".Last", envir=.GlobalEnv))
      {
        return(fun)
      }
    else
      {
        Last <- get(".Last", envir=.GlobalEnv)
        newfun <- function(...)
          {
            fun()
            Last()
          }
        return(newfun)
      }
  }

