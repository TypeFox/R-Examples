plot.safidesign <- function(x, runs = NULL, ...) {
  s.d <- x
  if (is.null(runs)) runs <- 1:nrow(s.d$DoE)
  if (class(s.d) != "safidesign") 
    stop("object class is not safidesign")
  l <- s.d$d.f
  m <- accessSafiDesign(s.d, n.timepoints = rep(1000, l))
  for (i in 1:length(m))
    m[[i]] <- m[[i]][runs,,drop=FALSE]
  if (l > 1) {
    nc <- round(sqrt(l))
    nl <- ceiling(l/nc)
    mfrow <- c(nc, nl)
    par(mfrow = mfrow)
  }
  for (i in 1:l) {
    matplot(t(m[[i]]), type = "l", main = s.d$variable.names[i], xaxt = "n", xlab = "input function domain", 
            ylab = "functional input", ylim = c(-1,1), ...)
    axis(1, at = c(0, 1000), labels = c(0, 1))
  }
  par(mfrow = c(1, 1))
}