`print.gpava` <-
function(x, ...)
{
  #x... object of class pava
  cat("\nCall:\n")
  print(x$call)
  cat("\nPAVA results:\n")
  res.table <- t(data.frame(x$z, round(x$x,2)))
  rownames(res.table) <- c("Predictor", "Fitted Values")
  colnames(res.table) <- 1:length(x$z)
  print(t(res.table))
  cat("\n")
  invisible(x)
}

