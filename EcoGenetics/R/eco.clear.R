#' Clearing the working environment, maintaining only the specified objects
#' 
#' @description This function removes all the elements of the working environment,
#' with the exception of those included in the argument of the function. 
#' Hidden elements can also be removed by setting all = TRUE.
#' @param ... Objects to retain.
#' @param all Remove also hidden elements? Default FALSE.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' ls()
#' eco.clear(eco)
#' ls()
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export

eco.clear <- function(..., all = FALSE) {
  
  clean.names <- as.character(match.call())[-1]
  
  env <- parent.frame()
  cuales <- ls(envir = env, all.names = all) 
  cuales  <- cuales %in% clean.names
  rm(list = ls(envir = env, all.names = all)[!cuales], envir = env)
  
}
