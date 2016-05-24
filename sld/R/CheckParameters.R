sl.check.pars <- function(pars) {
  if (length(pars) != 3) {
    ret <- FALSE
    warning("pars needs 3 elements, alpha, beta, delta\n")
  } else {
    ret <- TRUE # then break this is needed
    if (pars[2] <= 0) {
      ret <- FALSE
      warning("Negative or zero beta\n")
    } 
    if ((pars[3] < 0)|(pars[3] > 1)) {
      ret <- FALSE 
      warning("delta (skewing parameter) needs to be in the range [0,1]\n")
    } 
	}	
ret
}
