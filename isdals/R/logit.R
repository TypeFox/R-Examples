logit <- function(p) {
  if (!is.numeric(p))
    stop("input must be numeric")
  
  return(log(p) - log(1-p))
}
