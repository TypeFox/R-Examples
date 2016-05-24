#' @export
coef.rem.fi <- function(object,...){
  vcov <- solvecov(object$hessian)$inv
  coeff <- data.frame(estimate=coef(object$mr),
                      se=sqrt(diag(vcov)))
  return(coeff)
}
