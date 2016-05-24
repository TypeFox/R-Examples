print.t3way <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  df <- data.frame(value = c(x$Qa, x$Qb, x$Qc, x$Qab, x$Qac, x$Qbc, x$Qabc), 
                   p.value = c(x$A.p.value, x$B.p.value, x$C.p.value, x$AB.p.value, x$AC.p.value, x$BC.p.value, x$ABC.p.value))
  rownames(df) <- c(x$varnames[2], x$varnames[3], x$varnames[4], 
                    paste0(x$varnames[2], ":", x$varnames[3]),paste0(x$varnames[2], ":", x$varnames[4]),paste0(x$varnames[3], ":", x$varnames[4]),
                    paste0(x$varnames[2], ":", x$varnames[3], ":", x$varnames[4]))
  print(df)
  cat("\n")
}
