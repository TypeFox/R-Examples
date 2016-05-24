"pwm.beta2alpha" <-
function(pwm) {
  nmom <- length(pwm)
  nm1  <- nmom - 1
  otherpwm <- vector(mode="numeric", length=nmom)
  for(r in 0:nm1) {
     otherpwm[r+1] <- sum(sapply(0:r, function(k) { return((-1)^k*choose(r,k)*pwm[k+1]) }))
  }
  names(otherpwm) <- sapply(0:nm1, function(k) {
                            return(gsub("$", k, "Alpha", perl=TRUE)) } )
  return(otherpwm)
}


"pwm.alpha2beta" <-
function(pwm) {
  beta <- pwm.beta2alpha(pwm)
  names(beta) <- sapply(0:(length(beta)-1), function(k) {
                        return(gsub("$", k, "Beta", perl=TRUE)) } )
  return(beta)
}


