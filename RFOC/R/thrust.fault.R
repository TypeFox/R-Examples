`thrust.fault` <-
function(anim= seq(from=0, to=1, by=.1), KAPPA = 2,  Light=c(45,45))
  {
    if(missing(anim)) {  anim= seq(from=0, to=1, by=.1)  }
    if(missing(KAPPA)) {    KAPPA = 2 }
    if(missing(Light)) { Light=c(45, 45) }

    normal.fault(135, anim=anim, KAPPA=KAPPA, Light=Light)

  }

