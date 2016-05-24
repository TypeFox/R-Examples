# see coef.ds for documentation
#' @export
coef.io.fi <-function(object,...){
  if(length(grep("gam",as.character(object$model)))==0){
    vcov <- solvecov(object$hessian)$inv
  }else{
    vcov <- object$hessian
  }

  coeff <- data.frame(estimate=coef(object$mr),
                      se=sqrt(diag(vcov)))

  return(coeff)
}
