dataTrans <- function(y,mu,direction="backward",geoMean=NULL) {
  # Box-Cox transformation
  if (tolower(substr(direction,1,1))=="f") {
    if (any(y<=0)) stop("Observations must be strictly positive to be Box-Cox transformed")
    if (is.null(geoMean)) geoMean <- exp(mean(log(y)))
    if (mu==0) {
      return(geoMean*log(y))
    } else {
      return((y^mu-1)/mu/(geoMean^(mu-1)))
    }
  }
  # Inverse Box-Cox transformation 
  if (tolower(substr(direction,1,1))=="b") {
    if (is.null(geoMean)) stop("Geometric mean must be given for backward Box-Cox transformation")
    if (mu==0) {
      return(exp(y/geoMean))
    } else {
      return((1+y*mu*(geoMean^(mu-1)))^(1/mu))
    }
  }
}
