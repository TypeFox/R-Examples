groot <- function(S2, D2bool) {
  if(is.na(S2)) return(NA) ## can occur in g CI computation from condS2 CI computation
  ## ELSE
  if (D2bool) {
    if (S2<0.5) {
      ## message.redef("(!) Impossible S2 value in groot(S2)")
      return(NaN) ## not possible, _2_S2>=1, see Migraine doc
    }
    S2threshold <- (143 + 19*sqrt(57))/16
    if (S2<S2threshold) {
      machin <- (S2^2*(63 - 2*S2) + 3*sqrt(3)*sqrt(S2^3*(4 + (143 - 8*S2)*S2)))
      machin <- machin^(1/3)
      return((8*S2 - (2^(4/3)*(-3 + S2)*S2)/machin - 2^(2/3)*machin)/(6*S2))
    }
    ## else S2 > threshold
    theta <- acos(((63 - 2*S2)*S2^2)/(2*sqrt(((-3 + S2)*S2)^3)))
    return(4/3 - (2*((-3 + S2)*S2)^(1/2)*cos(theta/3))/(3*S2))
  } else {
    return((1+2*S2-sqrt(1+8*S2))/(2*S2))
  }
}
