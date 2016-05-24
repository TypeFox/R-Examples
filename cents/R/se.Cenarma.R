#se.Cenarma.R
se.Cenarma <- function(obj, nBoot=250){
  phi <- numeric(nBoot)
  p <- obj$p
  for (iBoot in 1:nBoot){
    y <- boot.Cenarma(obj)$y
    outi <- try(cenarma(y, p=p), silent=TRUE)
    if (identical(class(outi), "try-error")) phi[iBoot] <- NA
    else phi[iBoot] <- coef(outi$outarima)["ar1"]
  }
  sd(phi, na.rm=TRUE)
}
