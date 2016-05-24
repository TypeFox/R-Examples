"asHSdeCani" <- function( alpha, tk, oneOrTwoSided, gamma) {
  if ( gamma == 0 ) {
    return (alpha * tk)
  } else {
    return ((alpha/oneOrTwoSided) * ((1 - exp(-gamma * tk)) / (1 - exp(-gamma))))
  }
}


"asOBF" <- function( alpha, tk, oneOrTwoSided ) {
  # additonalParameters are not used here
  2 * (1 - pnorm ((qnorm(1 - (alpha / oneOrTwoSided)/2)) / sqrt(tk)))
}


"asPocock" <- function( alpha, tk, oneOrTwoSided ) {
  # additonalParameters are not used here
  (alpha/oneOrTwoSided) * log(1 + ( exp(1) - 1 ) * tk)
}


"asPowerFamily" <- function( alpha, tk, oneOrTwoSided, delta ) {
  (alpha/oneOrTwoSided)*(tk^delta)
}
