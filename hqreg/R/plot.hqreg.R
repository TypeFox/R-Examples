plot.hqreg <- function(x, xvar = c("lambda", "norm"), log.x =FALSE, nvars = TRUE, alpha = 1, ...)
{
  xvar <- match.arg(xvar)
  if (xvar == "norm" && log.x == TRUE) stop("set log.x to be FALSE when plotting against the L1 norm")
  nonzero <- which(rowSums(abs(x$beta))!=0)
  Y <- x$beta[nonzero[-1],]
  p <- nrow(Y)
  if (xvar == "lambda") {
    X <- x$lambda
    if (log.x) X <- log(X)
  } else {
    X <- colSums(abs(Y))
  }
  
  if (xvar == "lambda") {
    xlab <- if (log.x) expression(log(lambda)) else expression(lambda)
    xlim <- rev(range(X))
  } else {
    xlab <- expression(group("||", hat(beta), "||")[1])
    xlim <- range(X)
  }
  plot.args <- list(x=X, y=seq(X), ylim=range(Y), xlab=xlab, ylab="", 
                    type="n", xlim=xlim)
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  if (!is.element("ylab", names(new.args))) mtext(expression(hat(beta)), side=2, cex=par("cex"), line=3, las=1)
  
  cols <- hcl(h=seq(15, 375, len=max(4, p+1)), l=60, c=150, alpha=alpha)
  cols <- if (p==2) cols[c(1,3)] else cols[1:p]  
  line.args <- list(col=cols, lwd=1+2*exp(-p/20), lty=1)
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- X
  line.args$y <- t(Y)
  do.call("matlines",line.args)
  abline(h=0)
  if (nvars) {
    n.s <- apply(coef(x, lambda=x$lambda)!=0, 2, sum)-1
    axis(3, at=X, labels=n.s, tick=FALSE, line=-0.5)
  }
}

