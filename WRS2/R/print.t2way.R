print.t2way <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  if (!is.na(x$Qa)) {
    df <- data.frame(value = c(x$Qa, x$Qb, x$Qab), p.value = c(x$A.p.value, x$B.p.value, x$AB.p.value))
  } else {
    df <- data.frame(p.value = c(x$A.p.value, x$B.p.value, x$AB.p.value))
  }
  rownames(df) <- c(x$varnames[2], x$varnames[3], paste0(x$varnames[2], ":", x$varnames[3]))
  print(round(df,4))
  cat("\n")
}
