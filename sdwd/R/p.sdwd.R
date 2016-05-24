predict.sdwd = function(object, newx, s=NULL, 
    type=c("class", "link"), ...) {
  type = match.arg(type)
  b0 = t(as.matrix(object$b0))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[ , lamlist$left, drop=FALSE] %*% 
            Diagonal(x=lamlist$frac) +
            nbeta[ , lamlist$right, drop=FALSE] %*% 
            Diagonal(x=1-lamlist$frac)
    dimnames(nbeta) = list(vnames, paste(seq(along=s)))
  }
  nfit = as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  switch(type, link=nfit, class=ifelse(nfit > 0, 1, -1))
} 

print.sdwd = function(x, digits=max(3, 
    getOption("digits")-3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df=x$df, Lambda=signif(x$lambda, digits)))
} 


plot.sdwd = function(x, xvar=c("norm", "lambda"), 
    color=FALSE, label=FALSE, ...) {
  beta = x$beta
  lambda = x$lambda
  df = x$df
  xvar = match.arg(xvar)
  ##beta should be in 'dgCMatrix' format
  which = nonzero(beta)
  beta = as.matrix(beta[which, ])
  xvar = match.arg(xvar)
  switch(xvar, 
    norm = {
      index = colSums(abs(beta))
      iname = "L1 Norm"
           }, 
    lambda = {
      index = log(lambda)
      iname = "Log Lambda"
           }
  )
  xlab = iname
  ylab = "Coefficients"
  dotlist = list(...)
  type = dotlist$type
  if (is.null(type)) {
    if (color == FALSE) 
      matplot(index, t(beta), lty=1, xlab=xlab, ylab=ylab, 
              type="l", pch=500, 
              col=gray.colors(12, start=0.05, 
                  end=0.7, gamma=2.2), ...) else matplot(index, 
      t(beta), lty=1, xlab=xlab, ylab=ylab, 
      type="l", pch=500, ...)
  } else matplot(index, t(beta), lty=1, xlab=xlab, ylab=ylab, ...)
    atdf = pretty(index)
  prettydf = trunc(approx(x=index, y=df, xout=atdf, rule=2)$y)
  axis(3, at=atdf, labels=prettydf, cex.axis=1, tcl=NA)
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex=0.5, pos=pos)
  }
} 
