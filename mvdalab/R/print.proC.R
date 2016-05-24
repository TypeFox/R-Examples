print.proC <- function(x, ...) {
  ans <- rbind(x$Anova.MSS, x$Anova.ESS, x$Anova.TSS)
  dimnames(ans) <- list(NULL, "Sum of Squares")
  ans <- data.frame(Source = c("Model", "Error", "Total"), ans)
  print(list(`Rotation Matrix` = x$Rotation.Matrix, `Analysis of Variance` = ans))
}

