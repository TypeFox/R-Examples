simCorVars <- function(n, r, plot = TRUE){
  if(!is.integer(n)) n <- as.integer(n)
  if(length(n) > 1){
    n <- n[1]
    warning("length(n) > 1: only the first entry is used.")
  }
  stopifnot(is.numeric(r))
  stopifnot(abs(r) <= 1)
  if(length(r) > 1){
    r <- r[1]
    warning("length(r) > 1: only the first entry is used.")
  }
  x1 <- rnorm(n)
  y1 <- rnorm(n)
  a <- sqrt((1-r)/(1+r))
  if(abs(r) == 1){
    x <- x1
    y <- sign(r)*x1
  }else{
    x <- x1 + a*y1
    y <- x1 - a*y1
  }
  if(plot) plot(x, y, main = paste("Correlation =", r),
                xlab = "Variable 1", ylab = "Variable 2")
  data.frame(Var1 = x, Var2 = y)
}
