# ExtremeBounds - Extreme Bounds Analysis in R
# Author: Marek Hlavac

eba <-
function(formula = NULL, data, y = NULL, free=NULL, doubtful=NULL, focus=NULL, k=0:3, mu=0, level=0.95, vif=NULL, 
         exclusive=NULL, draws=NULL, reg.fun=lm, se.fun=NULL, include.fun=NULL, weights=NULL, ...) {
  
  cl <- match.call()
    
  return(.eba.wrap(formula=formula, data=data, y=y, free=free, doubtful=doubtful, focus=focus, 
                     k=k, mu=mu, level=level, vif=vif, exclusive=exclusive, draws=draws,
                     reg.fun=reg.fun, se.fun=se.fun, include.fun=include.fun, weights=weights, cl = cl, ...))
}

hist.eba <- function(x, variables=NULL, col="gray", freq=FALSE, main=NULL,
                     mu.show=TRUE, mu.col="red", mu.lwd=2, mu.visible=TRUE,
                     density.show=TRUE, density.col="blue", density.lwd=2, density.args=NULL, 
                     normal.show=FALSE, normal.col="darkgreen", normal.lwd=2, normal.weighted=FALSE, 
                     xlim=NULL, ylim=NULL, ...) {
  
  cl <- match.call()
  
  return(.hist.eba.wrap(x=x, variables=variables, col=col, freq=freq, main=main,
                        mu.show=mu.show, 
                        mu.col=mu.col, mu.lwd=mu.lwd, mu.visible=mu.visible,
                        density.show=density.show,
                        density.col=density.col, density.lwd=density.lwd, density.args=density.args,
                        normal.show=normal.show, normal.col=normal.col, normal.lwd=normal.lwd,
                        normal.weighted=normal.weighted, xlim=xlim, ylim=ylim, cl=cl, ...))

}

print.eba <- function(x, digits = 3, ...) {
  return(.print.eba.wrap(x=x, digits=digits, ...)) 
}

summary.eba <- function(object, ...) {
  return(object)
}

coefficients.eba <- function(object, ...) {
  return(object$coefficients)
}
