## plot method for unidimensional scaling

plot.uniscale <- function(x, main, pch = 19, ...) {
  
  if (missing(main)) main <- "Configurations Unidimensional Scaling" else main <- main
  plot(x$conf, rep(0, length(x$conf)), axes = FALSE, ann = FALSE, pch = pch, type = "o", ylim = c(-0.2, 0.8), ...)
  title(main)
  text(x$conf, rep(0, length(x$conf)) + 0.05, names(x$conf), srt = 90, adj = c(0,0.5))
}

