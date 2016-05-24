regrresidplot <- function(x, y, resid.plot=FALSE, fit.line=TRUE,
                          lm.object=lm(y~x),
                          x.name=names(lm.object$model)[2],
                          col=trellis.par.get()$plot.symbol$col,
                          col.yhat=NULL,
                          col.fit="gray80",
                          col.resid="gray40",
                          ...) {
  xyplot(y ~ x, resid.plot=resid.plot, fit.line=fit.line,
         lm.object=lm.object, x.name=x.name,
         col.fit=col.fit, col.resid=col.resid, ...,
         panel=function(x, y, ..., resid.plot, fit.line, lm.object, x.name) {
           panel.xyplot(x, y, ..., col=col, pch=19)
           if (fit.line) {
             newdata <- data.frame(x=cplx(101))
             names(newdata)[1] <- x.name
             newdata$y <- predict(lm.object, newdata=newdata)
             panel.lines(newdata[[1]], newdata$y, col=col.fit, ...)
           }
           yhat <- fitted(lm.object)
           if (!is.null(col.yhat))
             panel.points(x, yhat, col=col.yhat, ..., pch=20)
           panel.residSquare(x, y, yhat=yhat,
                             resid.plot=resid.plot, col=col.resid, ...)
         })
}

panel.residSquare <- function(x, y, yhat, resid.plot=FALSE, col="black", ...) {
  if (resid.plot==FALSE) return()
  if (resid.plot=="square") {
    xright <- x + convertUnit(unit(abs(yhat-y), "native"),
                              axisFrom="y",  typeFrom="dimension",
                              unitTo="native",
                              axisTo="x",    typeTo="dimension",
                              valueOnly=TRUE)
    panel.rect(xleft=x, ybottom=y, xright, ytop=yhat, border=col, ...)
  }
  else
    panel.segments(x, y, x, yhat, col=col, ...)
}

cplx <- function(length) {
  xlim <- current.panel.limits()$xlim
  scale(1:length, 1, (length-1)/diff(xlim)) + xlim[1]
}
