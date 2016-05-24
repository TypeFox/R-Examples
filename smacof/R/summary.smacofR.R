`summary.smacofR` <-
function(object, ...)
{
  cat("\n")
  cat("Subjects configurations (rows):\n")
  print(round(object$conf.row,4))
  cat("\n")
  cat("Objects configurations (columns):\n")
  print(round(object$conf.col,4))

  cat("\n\n")
  cat("Stress per point rows:\n")
  spp.perc.row <- object$spp.row/sum(object$spp.row)*100
  sppmat.row <- cbind(sort(object$spp.row), sort(spp.perc.row))
  colnames(sppmat.row) <- c("SPP","SPP(%)")
  print(round(sppmat.row, 4))
  cat("\n")

  cat("Stress per point columns:\n")
  spp.perc.col <- object$spp.col/sum(object$spp.col)*100
  sppmat.col <- cbind(sort(object$spp.col), sort(spp.perc.col))
  colnames(sppmat.row) <- c("SPP","SPP(%)")
  print(round(sppmat.row, 4))
  cat("\n") 
}



