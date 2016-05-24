#' Retrieves available symbols in the specified environment.
#' 
#' @param mode
#'    The mode of misspelled name.
#' @param envir
#'    The base environment to search variables.
#' @importFrom utils installed.packages
getNames <- function(mode, envir=.GlobalEnv){
   ## library not found
   if (mode == "lib") {
      packages <- installed.packages()
      return(attr(packages, "dimnames")[[1L]])
   }
   
   ## not an exported object
   if (mode == "export") {
      return(ls(envir=envir))
   }
   
   current <- envir
   table <- character(0L)
   while (!identical(current, emptyenv())) {
      variables <- ls(envir=current)
      variables <- variables[isVariableName(variables)]
      table <- union(table, variables)
      current <- parent.env(current)
   }
   
   ## FIXME: Bad performance
#    if (mode == "fun") {
#       table <- table[sapply(table, exists, envir=envir, mode="function")]
#    }
   
   sort(table)
}

#' Checks if the given \code{name} is valid as a variable name for R.
#' 
#' @param name
#'    A character vector to check.
isVariableName <- function(name) {
   make.names(name) == name
}
