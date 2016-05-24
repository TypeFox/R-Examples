#' Object Summary
#' 
#' Generic \code{summary} function for the \code{aovrasch} object
#'  
#' @param object      \code{aovrasch} object
#' @param ...         Additional arguments affecting the summary produced.
#' 
#' @author 
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#' Jan Steinfeld \email{jan.steinfeld@@univie.ac.at}
#' 
#' @importFrom stats printCoefmat
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # simulate Rasch model based data
#' # 100 persons, 20 items,
#' dat <- simul.rasch(100, items = seq(-3, 3, length.out = 20))
#' # reshape simulated data into 'long' format with balanced assignment 
#' # of examinees into two subgroups.
#' dat.long <- reshape.rasch(dat, group = rep(0:1, each = nrow(dat) / 2))
#' # apply three-way analysis of variance with mixed classification for testing the Rasch model.
#' res <- aov.rasch(dat.long) 
#' summary(res)
#' }
summary.aovrasch <- function(object, ...) {
  
  cat("\nThree-way analysis of variance with mixed classification \n\n")
  
  printCoefmat(object[4, , drop = FALSE], has.Pvalue = TRUE, cs.ind = 0)
  
  cat("\n")
  
  # Warning: Statistically significant Main Effect A
  if (object[1, "Pr(>F)"] < .05) {
    
    warning(paste0("Main effect A (group) is statistically significant, ", 
                   "F(1, ", object[2, "DF"], ") = ", formatC(object[1, "F value"], format = "f", digits = 3),  
                   ", p = ", formatC(object[1, "Pr(>F)"], format = "f", digits = 3),  ", i.e. results may not be trustworthy."))
    
  }
  
}