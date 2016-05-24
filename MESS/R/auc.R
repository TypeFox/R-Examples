auc <-
function(x, y, from = min(x), to = max(x), type=c("linear", "spline"), ...) 
{
    type <- match.arg(type)
  
    if (length(x) != length(y)) 
        stop("x and y must have the same length")
    if (length(unique(x)) < 2) 
        return(NA)

    if (type=="linear") {    
      values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
    } else {
      res <- integrate(splinefun(x, y, method="natural"), lower=from, upper=to)$value
    }
    res
}


