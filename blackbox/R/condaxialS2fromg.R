condaxialS2fromg <- function(gv, D2bool) {
  if (D2bool) {
    return((1+gv)/((2-gv)*(1-gv)^2))
  } else return((1+gv)/((1-gv)^2))
}
