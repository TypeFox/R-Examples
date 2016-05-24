setClass(
  Class = "T.cont",
  contains = "T.value"
)


setMethod(
  f = "panel",
  signature = "T.cont",
  definition = function(object, x, y) {
    ## call panel for T.value object
    callNextMethod(object, x, y)
    if(object@data$storeData) {
      dftab <- object@data$dftab
      coordsx <- object@data$coordsx
      coordsy <- object@data$coordsy
    } else {
      dftab <- eval(object@data$dftab, envir = sys.frame(object@data$frame))
      coordsx <- eval(object@data$coordsx, envir = sys.frame(object@data$frame))
      coordsy <- eval(object@data$coordsy, envir = sys.frame(object@data$frame))
    }
    
    dftab <- dftab / sum(dftab)
    
    f1 <- function(x, w) {
      w1 <- weighted.mean(w, x)
      w <- (w - w1)^2
      w2 <- sqrt(weighted.mean(w, x))
      return(c(w1, w2))
    }
    
    if(object@g.args$meanX) {
      w <- t(apply(dftab, 2, f1, w = coordsy))
      panel.points(x = coordsx, y = w[, 1], pch = 20, cex = 1.5, col = "black")
      panel.segments(coordsx, w[, 1] - w[, 2] , coordsx, w[, 1] + w[, 2], col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty, lwd = object@adeg.par$plines$lwd)
    }
    
    if(object@g.args$meanY) {
      w <- t(apply(dftab, 1, f1, w = coordsx))
      panel.points(x = w[, 1], coordsy, pch = 20, cex = 1.5, col = "black")
      panel.segments(w[, 1] - w[, 2], coordsy, w[, 1] + w[, 2], coordsy, col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty, lwd = object@adeg.par$plines$lwd)
    }
    
    coordsx <- coordsx[col(as.matrix(dftab))]
    coordsy <- coordsy[row(as.matrix(dftab))]
    if(object@g.args$ablineX)
      panel.abline(reg = lm(coordsy ~ coordsx, weights = as.vector(as.matrix(dftab))), col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty, lwd = object@adeg.par$plines$lwd)
    
    if(object@g.args$ablineY) {
      w <- coefficients(lm(coordsx ~ coordsy, weights = as.vector(as.matrix(dftab))))
      if(w[2] == 0)
        panel.abline(h = w[1], col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty, lwd = object@adeg.par$plines$lwd)
      else
        panel.abline(c(-w[1] / w[2], 1 / w[2]), col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty, lwd = object@adeg.par$plines$lwd)
    } 
  })
