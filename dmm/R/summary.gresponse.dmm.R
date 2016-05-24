summary.gresponse.dmm <-
function(object, ...)
# summary.gresponse.dmm() - summarize a gresponse.dmm object (less brief output)
{
  cat("Call:\n")
  print(object$call)
  cat("\nPredicted response to selection using component(s) ",object$effects,"\n")

  cat("\nGenetic selection differentials achieved by given psd:\n\n")
  cat("Overall:\n")
  print(object$overall,digits=object$digits)
  cat("\n")
  cat("Sex specific:\n")
  print(object$sex$gsd,digits=object$digits)
  print(object$sex$gsdsum,digits=object$digits)
  print(object$sex$psd,digits=object$digits)
  cat("\n")
  cat("Path specific:\n")
  print(object$path$gsd,digits=object$digits)
  print(object$path$gsdsum,digits=object$digits)
  print(object$path$psd,digits=object$digits)


# cat("\nGenetic selection differentials achieved by a unit psd on each trait:\n\n")
# print(object$ugsd,digits=object$digits)
# cat("\n")

# cat("Directional selection gradient:\n\n")
# print(object$dsg,digits=object$digits)
# cat("\n")
#
# cat("Genetic var/covariance matrix:\n\n")
# print(object$gcov,digits=object$digits)
# cat("\n")
#
# cat("Phenotypic var/covariance matrix:\n\n")
# print(object$pcov,digits=object$digits)
# cat("\n")
# 
# cat("Phenotypic selection differentials:\n\n")
# print(object$psd,digits=object$digits)
# cat("\n")
}
