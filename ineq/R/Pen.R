Pen <- function(x, n = rep(1, length(x)), group = NULL,
  scaled = TRUE, abline = TRUE, add = FALSE, segments = NULL,  
  main = "Pen's Parade", ylab = NULL, xlab = NULL, 
  col = NULL, lwd = NULL, las = 1, fill = NULL, ...) 
{
  ina <- !is.na(x)    
  n <- n[ina]
  x <- as.numeric(x)[ina]
  o <- order(x)
  x <- x[o]
  n <- n[o]

  if(!is.null(group)) {
    if(!is.factor(group)) group <- as.factor(group)
    if(is.null(segments)) segments <- TRUE
    if(is.null(fill) && segments) {
      fill <- grey(seq(0.3 ^ 2.2, 0.9 ^ 2.2, length = length(levels(group))) ^ (1/2.2))
      fill <- fill[as.numeric(group)]
    }
  } else {
    if(is.null(segments)) segments <- FALSE
  }
  
  if(is.null(col)) {
    if(segments | !is.null(fill)) col <- 1
      else col <- 4
  }
  if(is.null(lwd)) {
    if(segments | !is.null(fill)) lwd <- 1
      else lwd <- 2
  }
  
  if(scaled) x <- x/mean(x)

  if(is.null(ylab)) {
    if(scaled) ylab <- expression(x[(i)]/bar(x))
      else ylab <- expression(x[(i)])
  }

  if(is.null(xlab)) {
    if(identical(all.equal(n, rep(1, length(x))), TRUE)) xlab <- expression(i/n)
      else xlab <- expression(Sigma[i](n[(i)]/n))
  }
  
  n <- cumsum(c(0, n))/sum(n)   

  if(!add) plot(c(0, 1), c(0, max(x)), type = "n", main = main, ylab = ylab, xlab = xlab, 
      xaxs = "i", yaxs = "i", las = las, ...)

  if(segments) {
    lx <- length(x)
    col <- rep(col, length.out = lx)[o]
    fill <- rep(fill, length.out = lx)[o]
    for(i in 1:lx)
      polygon(c(n[i],n[i+1],n[i+1],n[i]), c(0, 0, x[i],x[i]), col = fill[i], border = col[i], lwd = lwd)
  } else {
    ln <- length(n)
    n2 <- c(rep(n[-ln], rep(2, ln-1)), n[ln])
    x2 <- c(0, rep(x, rep(2, ln-1)))
    lines(n2, x2, col = col, lwd = lwd)
    if(!is.null(fill)) polygon(c(n2, 1, 0), c(x2, 0, 0), col = fill, border = col)
  }
  
  if(abline & !add) abline(h = mean(x), lty = 3)
  box()
}
