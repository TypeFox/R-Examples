plot.mvdareg <- function (x, plottype = c("PE", "scoresplot", "loadingsplot", 
                             "loadingsplot2D", "T2", "Xresids", "coefsplot", "ap.plot", 
                             "weightsplot", "weightsplot2D", "acfplot"), ...) {
  object <- x
  plottype <- match.arg(plottype)
  plotFunc <- switch(plottype, PE = PE, scoresplot = scoresplot, loadingsplot = loadingsplot, 
                     loadingsplot2D = loadingsplot2D, T2 = T2, Xresids = Xresids, 
                     coefsplot = coefsplot, ap.plot = ap.plot, 
                     weightsplot = weightsplot, weightsplot2D = weightsplot2D, 
                     acfplot = acfplot)
  plotFunc(object, ...)
}

