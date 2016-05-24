Kinhomhat <-
function(X, r = NULL, ReferenceType = "", lambda = NULL, CheckArguments = TRUE) {

  if (CheckArguments)
    CheckdbmssArguments()
  
  # K intratype calls Kest with the best edge-effect correction and returns the values
  Kiintra <- function (X, r, lambda) {
    # Estimate intensity if it has not been provided
    if (is.null(lambda)) {
      lambda <- as.vector(density.ppp(X, sigma=bw.diggle(X), at="points"))
    }
    # Calculate Kinhom according to lambda
    return (Kinhom(X, r=r, correction="best", normpower=2, lambda=lambda))
  } 
  
  # K intra
  if (ReferenceType == "") {
    return (Kiintra(X, r, lambda))
  }
  # K intra for a single point type
  if (ReferenceType != "") {
    X.reduced <- X[X$marks$PointType==ReferenceType]
    return (Kiintra(X.reduced, r, lambda))
  }  

}
