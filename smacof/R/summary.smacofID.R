`summary.smacofID` <-
function(object, ...)
{
  cat("\n")
  cat("Group Stimulus Space (Joint Configurations):\n")
  print(round(object$gspace, 4))

  cat("\n\n")
  cat("Stress per point:\n")

  print(round(object$spp,2))
  cat("\n")

}


