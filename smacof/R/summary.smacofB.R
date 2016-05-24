`summary.smacofB` <-
function(object, ...)
{
  cat("\n")
  cat("Configurations:\n")
  print(round(object$conf,4))
  
  cat("\n\n")
  cat("Stress per point (in %):\n")

  #spp.perc <- object$spp/sum(object$spp)*100
  #sppmat <- cbind(sort(object$spp), sort(spp.perc))
  #colnames(sppmat) <- c("SPP","SPP(%)")
  print(round(object$spp,2))
  cat("\n")
}

