#' Create name for prcomp result
#' 
#' Create name for \code{\link[stats]{prcomp}} result
#' 
#' @param prcompResult output value from \code{\link[stats]{prcomp}} function
#' @param PC PC number 
#' @param explvar (logical) show explained variance (\%) or not
#' 
#' @return String
#' 
#' @examples
#' data(applejuice)
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' result <- prcomp(applejuice_uf) 
#' prcompname(result, 1)
#' 
#' @export
#' 
prcompname <- function(prcompResult, PC, explvar = TRUE) {
  # get information from prcompResult
  score <- prcompResult$x
  if (explvar == TRUE){
    varP <- round((prcompResult$sdev^2 / sum(prcompResult$sdev^2)) * 100, 
                  digits = 1)
    prcompname <- paste("PC ", PC, " (", varP[PC], "%)", sep = "")
  } else {
    prcompname <- paste("PC ", PC, sep = "")
  }
  return(prcompname)
}