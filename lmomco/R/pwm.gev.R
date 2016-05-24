"pwm.gev" <-
function(x,nmom=5,sort=TRUE) {
  z <- pwm.pp(x,A=-0.35,B=0,nmom=nmom,sort=sort)
  z$source <- "pwm.gev"
  return(z)
}
