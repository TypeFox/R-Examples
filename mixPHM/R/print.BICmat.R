`print.BICmat` <-
function(x,...)
# object of class BICmat
{
  dimnames(x$BICmat) <- list(x$method,paste("K =",x$K))
  cat("\n")
  cat("Survival distribution:",x$Sdist,"\n")
  cat("\n")
  print(x$BICmat)
  cat("\n")
}

