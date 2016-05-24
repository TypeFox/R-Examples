"Bpwm2ApwmRC" <-
function(Bpwm,para) {
  nmom <- NULL
  Abetas <- vector(mode="numeric")
  Bbetas <- vector(mode="numeric")
  
  if(is.list(Bpwm)) {
    Bbetas <- Bpwm$betas
  }
  else {
    Bbetas <- Bpwm
  }

  nmom <- length(Bbetas)
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
    Abetas[r-1+1] <- (r*Bbetas[r-1+1] - (1-zeta^r)*x.of.zeta)/(zeta^r*r)
  }
  z <- list(betas=Abetas,source="Bpwm2ApwmRC",XofZeta=x.of.zeta)
  return(z)
}

