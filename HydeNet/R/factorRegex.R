#' @name factorRegex
#' @importFrom stats terms
#' 
#' @title Produce Regular Expressions for Extracting Factor Names and Levels
#' @description A utility function to produce a regular expression that can
#'   separate factor names and factor levels in the \code{broom::tidy()$term}
#'   output.  At some point, this may have to become a method to accomodate
#'   different model types, but I haven't run into that problem yet.
#'   
#' @param fit a model object
#' 
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @examples 
#' data(PE, package = "HydeNet")
#' g6 <- glm(treat ~ d.dimer + angio, data=PE, family="binomial")
#' HydeNet:::factorRegex(g6)
#' 

factorRegex <- function(fit){
  fctr <- attributes(stats::terms(fit))$dataClasses
  if (any(fctr == "factor")){
    fctr <- names(fctr)[fctr == "factor"]
    fctr_regex <- paste0(fctr, collapse="|")
    fctr_regex <- gsub("[(]", "[(]", fctr_regex)
    fctr_regex <- gsub("[)]", "[)]", fctr_regex)
    fctr_regex <- gsub("[.]", "[.]", fctr_regex)
    fctr_regex <- paste0("(", fctr_regex, ")")
    return(fctr_regex)
  }
  else NULL
}




