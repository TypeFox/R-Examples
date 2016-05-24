"pairs.profile" <-
  function (x, labels = c(names(x), "Profile tau"), panel = lines, 
            invert = TRUE, plot.tau = TRUE, plot.trace = TRUE, plot.sketch = TRUE, 
            plot.ellipse = FALSE, level = 0.95, ...) 
{
  doaxis <- function(which, dolabel = TRUE) axis(which, labels = dolabel) # outer = TRUE, line = -0.5, labels = dolabel)
  setup <- function(x, y, ...) plot(range(x[!is.na(x)]), 
                                       range(y[!is.na(y)]), type = "n", axes = FALSE, ...)
  if (is.character(panel)) 
    panel <- get(panel, mode = "function")
  n <- length(x)
  if (plot.tau) 
    n <- n + 1
  oldpar <- par("oma", "mar", "cex", "tck", "mgp", "mex", 
                "mfrow")
  oldcex <- par("cex")
  CEX <- oldcex * max(7.7/(2 * n + 3), 0.6)
  par(mfrow = c(n, n), mgp = c(2, 0.8, 0), oma = rep(3, 4), 
      mar = rep(0.5, 4), tck = -0.03/n)
  on.exit({
    par(oldpar)
  })
  par(cex = CEX)
  if (length(labels) < n) 
    labels <- paste(deparse(substitute(x)), "[,", 1:n, "]", 
                    sep = "")
  if (par("pty") == "s") {
    dif <- diff(par("fin"))/2
    if (dif > 0) 
      par(omi = c(dif * n, 0, dif * n, 0) + par("omi"))
    else par(omi = c(0, (-dif) * n, 0, (-dif) * n) + par("omi"))
  }
  alltau <- unlist(lapply(x, function(x) x[[1]]), use.names = FALSE)
  order <- if (invert) 
    1:n
  else n:1
  for (i in order) {
    for (j in 1:n) {
      if (i <= length(x)) {
          icomp <- x[[i]]
	  ipars <- as.matrix(icomp[[2]])
      }
      if (j <= length(x)) {
          jcomp <- x[[j]]
	  jpars <- as.matrix(jcomp[[2]])
      }
      xx1 <- NA
      xx2 <- NA
      yy1 <- NA
      yy2 <- NA
      if (i <= length(x)) {
        yy1 <- ipars[, i]
        if (j <= length(x)) {
          xx1 <- ipars[, j]
          xx2 <- jpars[, j]
          yy2 <- jpars[, i]
        }
        else {
          xx1 <- icomp[[1]]
        }
      }
      else {
        yy1 <- jcomp[[1]]
        if (j <= length(x)) {
          xx1 <- jpars[, j]
        }
      }
      xx <- c(xx1, NA, xx2)
      yy <- c(yy1, NA, yy2)
      if (i <= length(x)) {
        if (j <= length(x)) 
          setup(xx, yy, ...)
        else setup(alltau, yy, ...)
      }
      else {
        if (j <= length(x)) 
          setup(xx, alltau, ...)
        else setup(alltau, alltau)
      }
      box()
      if (i == 1) 
        doaxis(3, j%%2 == 0)
      if (i == n) 
        doaxis(1, j%%2 == 1)
      if (j == 1) 
        doaxis(2, i%%2 == 0)
      if (j == n) 
        doaxis(4, i%%2 == 1)
      if (i != j) {
        if ((i <= length(x)) && (j <= length(x))) {
          if (plot.trace) 
            panel(xx, yy, ...)
          if (plot.sketch) 
            for (l in level) panel(ellipse(x, which = c(j, 
                                                i), level = l), ...)
          if (plot.ellipse && !is.null(fit <- attr(x, 
                                                   "original.fit"))) 
            for (l in level) panel(ellipse(fit, which = c(j, 
                                                  i), level = l), ...)
        }
        else if (plot.tau) 
          panel(xx, yy, ...)
      }
      else {
        par(usr = c(0, 1, 0, 1))
        text(0.5, 0.5, labels[i], cex = 1.5 * CEX)
      }
    }
  }
  invisible()
}
