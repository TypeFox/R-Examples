#' Function to get EC values given a model.
#'
#' @param x a model generated from the drm function of the drc package
#' @param response a numeric vector specifying at what response levels to generate an EC value. Defaults to 10, 50 and 90 (percent).
#' @param disp display the results?
#' 
#' @details This function is a wrapper for \code{\link[drc]{ED}}. 
#' @keywords internal
#'  
#' @return a data frame with the Estimate and standard error for each response
#'  value. 
#' @export
#' @author Zhian N. Kamvar
get_EC <- function(x, response = c(10, 50, 90), disp = TRUE){
  resnames <- c("Estimate", "SE")
  if (length(x) < 1){
    res <- matrix(as.numeric(NA), nrow = 1, ncol = 6)
  } else {
    res <- drc::ED(x, respLev = response, display = disp)
    res <- matrix(t(res), nrow = 1)
  }
  colnames(res) <- paste(resnames, rep(response, each = 2), sep = ".")
  return(data.frame(res))
}