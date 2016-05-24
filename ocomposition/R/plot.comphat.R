

plot.comphat <- function(x, ...){
  ci <- apply(x, 2, function(z) HPDinterval(as.mcmc(z), prob = 0.95))
  b <- barplot(apply(x, 2, mean), names.arg = 1:ncol(x), las = 1, ylim = c(0, max(c(ci))), col = "dark blue")
  segments(b, ci[1,], b, ci[2,], lwd = 5, col = "dark red")
  }