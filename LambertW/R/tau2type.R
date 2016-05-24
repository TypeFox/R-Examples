#' @rdname tau-utils
#' @description
#' \code{tau2type} guesses the type (\code{'s'}, \code{'h'}, \code{'hh'}) from the names
#' of \code{tau} vector; thus make sure \code{tau} is named correctly.
#' @export
#' 
#' 

tau2type <- function(tau) {
  
  if (!("gamma" %in% names(tau))) {
    tau["gamma"] <- 0
  }
  if (tau["gamma"] == 0) {
    if (any(grepl("_l|_r", names(tau)))) {
      type <- "hh"
    } else {
      type <- "h"
    }
  } else {
    type <- "s"
  }
  return(type)
}