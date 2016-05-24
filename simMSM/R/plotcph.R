plotcph <- function(m, ...){
  b <- m$coef
  ci <- confint(m)
  P <- length(b)
  plot(1:P, seq(max(ci), min(ci), length = P), xaxt = "n", xlab = "", 
       ylab = "Covariate effect size", type = "n", xlim = c(0.6, P + 0.4), 
       bty = "n", ...)
  abline(h = 0)
  axis(1, at = 1:P, names(b), las = 2)
  for(p in 1:P){
    x <- p + 0.3*c(-1, 1)
    y <- ci[p, ]
    polygon(x = c(x, rev(x)), y = c(y[1], y[1], y[2], y[2]), col = "lightgrey", border = NA)
    lines(x, rep(b[p], 2), lwd = 2)
  }
}