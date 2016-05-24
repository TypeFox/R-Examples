#' Computing Sampling error
#' 
#' @param sm sample mean
#' @param mu population mean
#' @export
#' @examples
#' sampling.error(sm=0.45, mu=0.5)
sampling.error <- function(sm,mu){
  return(sm-mu)
}

