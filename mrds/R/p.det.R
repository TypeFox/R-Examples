#' Double-platform detection probability
#'
#' Computes detection probability for detection function computed from
#' mark-recapture data with possibly different link functions.
#'
#' @param dpformula formula for detection function
#' @param dplink link function ("logit","loglog","cloglog")
#' @param dppars parameter vector
#' @param dpdata double platform data
#' @return vector of predicted detection probabilities
#' @author ?????
#' @importFrom stats model.frame model.offset
p.det <- function(dpformula, dplink, dppars, dpdata){
  fm <- dpformula
  modframe <- model.frame(fm, data=dpdata, drop.unused.levels=FALSE)
  dat <- model.matrix(fm, data=modframe)
  offsetval <- model.offset(modframe)

  if (length(dppars)>1){
    lpred <- dat %*% dppars
  }else{
    lpred <- dat * dppars
  }

  if(!is.null(offsetval)){
    lpred <- lpred+offsetval
  }

  p <- switch(dplink,
              loglog  = exp(-exp(lpred)),
              cloglog = 1-exp(-exp(lpred)),
              logit   = exp(lpred)/(1+exp(lpred))
     )
  return(as.vector(p))
}
