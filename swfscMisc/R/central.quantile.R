#' @title Central Quantile
#' @description Upper and lower values of central quantile
#' 
#' @param x numeric vector.
#' @param pct central percentile desired.
#' 
#' @return a two element vector giving the lower and upper quantiles.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' x <- runif(1000)
#' central.quantile(x)
#' central.quantile(x, pct = 0.75)
#' 
#' @importFrom stats quantile
#' @export
#' 
central.quantile <- function(x, pct = 0.95) {
  if((pct < 0) | (pct > 100)) return(NA)
  if(pct > 1) pct <- pct / 100
  lci <- (1 - pct) / 2
  uci <- 1 - lci
  cent.quant <- quantile(x, probs = c(lci, uci))
  names(cent.quant) <- c(paste("lower", pct, sep = "_"), paste("upper", pct, sep = "_"))
  cent.quant
}
