#' Deletes variables from dataset
#'
#' R implementation of the SPSS \code{DELETE VARIABLES} argument.
#' 
#' @usage xpssDeleteVariables(x, variables = NULL)
#' @param x input data.
#' @param variables atomic character or character vector with the name of the variables.
#' @return Output is the narrowed dataset.
#' @author Bastian Wiessner
#' @seealso Related Functions \code{\link{drop}} \code{\link{subset}} 
#' @examples 
#' data(fromXPSS)
#' xpssDeleteVariables(fromXPSS, variables = "V1")
#' @export
xpssDeleteVariables <- function(x, variables = NULL){
  #do meta check  
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  if(unique(!(is.element(variables,names(x))))){
    stop("variables arent in the dataset")
  }
  
  pos <- which(names(x) %in% variables)
  x <- x[c(-pos)]
  
  x <- applyAttributeDemerge(x)

  return(x)
}
