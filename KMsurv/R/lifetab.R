lifetab <- function (tis, ninit, nlost, nevent) {
  ## tis has length 1 more than other vectors
  
  Ypj <- c(ninit, ninit - cumsum(nlost + nevent)[-length(nevent)])
  Yj <- Ypj - nlost/2
  Sj <- cumprod(1 - nevent/Yj)

  qj <- nevent/Yj
  pj <- 1 - qj
  n <- length(Yj)
  Sj <- c(1, Sj[-n])
  fmj <- c(diff(-1 * Sj), NA)/diff(tis)
  hmj <- nevent / diff(tis) / (Yj - nevent/2)
  hmj[n] <- NA

  Sj.se <- c(0, Sj[-1] * sqrt(cumsum(nevent/Yj/(Yj - nevent))[-length(Sj)]))
  
  fmj.se <- Sj*qj/diff(tis) * sqrt(c(0,cumsum(qj/Yj/pj)[-n]) + (pj/Yj/qj))
  fmj.se[n] <- NA
  
  hmj.se <- sqrt(1 - (hmj * diff(tis)/2)^2) * sqrt(hmj^2/Yj/qj)
  hmj.se[n] <- NA
  
  data.frame(nsubs=Ypj, nlost=nlost, nrisk=Yj, nevent=nevent,
             surv=Sj, pdf=fmj, hazard=hmj,
             se.surv=Sj.se, se.pdf=fmj.se, se.hazard=hmj.se,
             row.names=paste(tis[-n-1], tis[-1], sep="-"))
}

