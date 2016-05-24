print.gof <- function(x, ..., digits = max(3, getOption("digits") - 3)){
  
  D <- x$D
  X2 <- x$X2
  df <- x$df
  
  cat("D  = ", round(D, digits = digits), ", df = ", df,
    ", P(>D) = ", pchisq(D, df, lower.tail = FALSE), '\n', sep = "")

  cat("X2 = ", round(X2, digits = digits), ", df = ", df,
    ", P(>X2) = ", pchisq(X2, df, lower.tail = FALSE), '\n', sep = "")
    
}
