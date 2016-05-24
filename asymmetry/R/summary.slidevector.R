summary.slidevector <-
function(object,...)
{
  cat("\n")
  cat("Configurations:\n")
  print(round(object$confi,4))

  cat("\n")
  cat("Slide-vector:\n")
  cat(round(object$slvec,4))

}
