seaRoll <- 
function(x, w = 10, plot = FALSE, pval = .05, mval = .5, 
         pchs = c(19, 21), xlab = NULL, ylab = NULL, ...) {
         
  # validate args
  if (!is.ts(x) || is.matrix(x))
    stop("'x' must be a vector of class 'ts'")
  if (w < 5)
    stop("window must be at least 5 years")
  
  # result for each window
  sr <- function(y, x1=x, w1=w, mval1=mval) {
    # set current window and get slope
    fr <- frequency(x1)
    x2 <- window(x1, start = y, end = c(y + w1 - 1, fr), extend = TRUE)
    sk <- seaKen(x2)
    
    # make sure enough data are present
    miss.ok <- sum(sk$miss >= 0.5) / fr < mval1
    N <- sum(!is.na(x2))
    if (N < 3 * fr || N < 10 || !miss.ok) {
      c(NA, NA, NA)
    } else {
      c(sk$sen.slope, sk$sen.slope.rel, sk$p.value)
    }
  }
  
  # combine windows
  sx <- start(x)[1]
  ex <- end(x)[1]
  ans <- t(sapply(sx:(ex - w + 1), function(i) sr(i)))
  rownames(ans) <- sx:(ex - w + 1)
  colnames(ans) <- c("sen.slope", "sen.slope.rel", "p.value")

  # plot if TRUE
  if (!plot) {
    ans
  } else {
    xlab <- if (is.null(xlab)) "" else xlab
    ylab <- if (is.null(ylab)) "" else ylab
    pch <- ifelse(ans[, "p.value"] < pval, pchs[1], pchs[2])
    plot(ans[, "sen.slope"] ~ rownames(ans), pch = pch, 
         xlab = xlab, ylab = ylab, ...)
  }
}
