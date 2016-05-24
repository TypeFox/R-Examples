#' @export 
plot.rankFD <- function (x, ...) {
  
  object <- x
  dots <- list(...)
  a <- object$plotting
  
  # default values
  args <- list(nf = a$nf, fac_names = a$fac_names, n.hypotheses = a$n.hypotheses,
               Descriptive.Factors = a$Descriptive.Factors, CI.method = a$CI.method)
  
  args[names(dots)] <- dots
  
  do.call(plotting, args = args)
}

#' @export
print.rankFD <- function (x, ...) {
  cat("\n", "Descriptive:", "\n", sep = "")
  print(x$Descriptive)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(x$Wald.Type.Statistic)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(x$ANOVA.Type.Statistic)
}
