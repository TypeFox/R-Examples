plot.ifile <- function(x, theta = seq(-4, 4, length = 100),
  subset = 1:nrow(x), ...) {

  if(is.character(subset))
    subset <- which(x$name == subset)
  b <- x$measure[subset]

  if(interactive()) {
    plot(rirf(b, theta), ...)
    readline("Enter to see next plot")
    plot(riif(b, theta), ...)
    readline("Enter to see next plot")
    plot(rief(b, theta), ...)
    readline("Enter to see next plot")
    plot(rtrf(b, theta), ...)
    readline("Enter to see next plot")
    plot(rtif(b, theta), ...)
    readline("Enter to see next plot")
    plot(rtef(b, theta), ...)
  }
}
