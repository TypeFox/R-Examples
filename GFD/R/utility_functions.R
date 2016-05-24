#' @export 
plot.GFD <- function (x, ...) {
  
  object <- x
  dots <- list(...)
  a <- object$plotting
  b <- object$Descriptive
  fac.names <- a$fac_names
  exist <- hasArg(factor) 
   
  if(length(fac.names) != 1){
    if(!exist){
      print("Please choose the factor you wish to plot (for interaction type something like group1:group2) and confirm by pressing 'Enter'")
      Faktor <- scan("", what="character")
      while(length(Faktor)==0){
        print("Please enter the name of the factor you wish to plot!")
        Faktor <- scan("", what="character")
      }
      } else {
       Faktor <- dots$factor
      } 
    x.label <- ""
  } else {
    Faktor <- fac.names
    x.label <- fac.names
  }
    
  match.arg(Faktor, fac.names)

  # default values
  args <- list(plot.object = a, descr.object = b, factor = Faktor,
               lwd =2, ylab = "Means", xlab = x.label, col = 1:length(fac.names), pch = 1:18, legendpos = "topright")
  
  args[names(dots)] <- dots
  
  do.call(plotting, args = args)
}

#' @export
print.GFD <- function(x, ...) {
  a <- x$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(x$WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(x$ATS)
}

#' @export
summary.GFD <- function (object, ...) {
  a <- object$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(object$WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(object$ATS)
}
