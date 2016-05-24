from2Ns2Tocanon <- function(parvec) { ## parvec is named c(2Nmu, latt2Ns2, 2Nm ou g)
  ## imp the [[]] remove the names, which are otherwise appended to the given names in c(twoNmu=th...)
  th <- parvec[["twoNmu"]]
  latt2Ns2 <- parvec[["latt2Ns2"]]
  if ("twoNm" %in% names(parvec)) {
    mig <- parvec[["twoNm"]]
    g <- groot(latt2Ns2/mig, D2bool=("2D" %in% blackbox.getOption("DemographicModel")) ) ## g=groot(sigma2|disp).
  } else if ("g" %in% names(parvec)) {
    g <- parvec[["g"]]
    mig <- latt2Ns2/condaxialS2fromg(g, D2bool=("2D" %in% blackbox.getOption("DemographicModel")))## 2Nm
  }
  return(c(twoNmu=th, twoNm=mig, g=g))
}
