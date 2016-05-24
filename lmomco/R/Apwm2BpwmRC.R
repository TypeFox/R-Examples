"Apwm2BpwmRC" <-
function(Apwm,para) {
  nmom <- NULL
  Abetas <- vector(mode="numeric")
  Bbetas <- vector(mode="numeric")
  
  if(is.list(Apwm)) {
    Abetas <- Apwm$betas
  }
  else {
    Abetas <- Apwm
  }

  nmom <- length(Abetas)
  zeta <- NULL
   
  if(is.na(para$zeta)) {
     warning("Zeta is missing in parameter object, suspect you do not have parameters for a right-censored distribution")
     return(NULL)
  }
  else {
     zeta <- para$zeta
  }
  if(zeta < 0 | zeta > 1) {
     warning("Bad zeta: 0 <= zeta <= 1")
     return(NULL)
  }

  x.of.zeta <- par2qua(zeta,para)

  for(r in seq(1,nmom)) {
    Bbetas[r-1+1] <- (zeta^r*r*Abetas[r-1+1] + (1-zeta^r)*x.of.zeta)/r
  }
  z <- list(betas=Bbetas,source="Apwm2BpwmRC",XofZeta=x.of.zeta)
  return(z)
}
