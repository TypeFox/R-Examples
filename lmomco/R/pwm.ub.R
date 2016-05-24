"pwm.ub" <-
function(x,nmom=5,sort=TRUE) {
  z <- pwm(x,nmom=nmom,sort=sort)
  z$source <- "pwm.ub"
  return(z)
}
