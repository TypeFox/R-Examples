summary.flim <-
function(object, ...) {
  cat("##############################################################")
  cat("\nCall:\n")
  print(object$call) 
  total.ind <- length(unique(object$df[, 1]))
  CountObserved <- function(vector) sum(vector==1)/total.ind
  obs.fractions <- tapply(object$df[, "obs.type"], object$df[, 2],
                          FUN = "CountObserved")
  ofrac <- as.data.frame(rbind(object$times, round(obs.fractions, 2)))
  rownames(ofrac) <- c("Timepoint", "Fraction")
  cat("\nObserved fractions on each time point:\n")
  print(round(obs.fractions, 2))
  cat("##############################################################\n")
  cat("\nRegression models for increments:\n")
  cat("Response[t+1] - Response[t]\n")
  cat("\n")
  print(flimList(object))
  cat("\n")
}
