#' @title print.damuth
#'  
#' @description S3 method for class \code{damuth}
#'
#' @details
#' See Examples
#' 
#' @param x an object of class \code{damuth}
#' @param ... arguments to be passed to methods
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(arth$spp, arth$count, arth$mass^0.75)
#' ebar1 <- ebar(esf1)
#' print(ebar1)
#' 
#' @return Returns the object silently
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow

print.damuth <- function(x,...) {
  cat('Abundance metabolic rate relationship ranging from \n')
  
  cat(sprintf('n: [%s, %s] \ne: [%s, %s] \n', 
              min(x[['n']], na.rm=TRUE), 
              max(x[['n']], na.rm=TRUE), 
              min(x[['e']], na.rm=TRUE), 
              max(x[['e']], na.rm=TRUE)))
  
  invisible(x)	
}