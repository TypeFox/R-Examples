plot.grpreg <- function(x, alpha=1, legend.loc, log.l=FALSE, norm=FALSE, ...) {
  if (norm) {
    Y <- predict(x, type="norm")
    if (any(x$group==0)) Y <- Y[-1,]
    nonzero <- which(apply(abs(Y), 1, sum)!=0)
    Y <- Y[nonzero,]
    g <- 1:nrow(Y)
  } else {
    beta <- if (length(dim(x$beta))==3) matrix(x$beta[,-1,,drop=FALSE], ncol=dim(x$beta)[3]) else x$beta[-1,,drop=FALSE]
    penalized <- which(x$group!=0)
    nonzero <- which(apply(abs(beta),1,sum)!=0)
    ind <- intersect(penalized, nonzero)
    Y <- beta[ind,,drop=FALSE]
    g <- as.numeric(as.factor(x$group[ind]))
  }
  p <- nrow(Y)
  l <- x$lambda
  n.g <- max(g)
  
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)
  
  plot.args <- list(x=l, y=1:length(l), ylim=range(Y), xlab=xlab, ylab="", type="n", xlim=rev(range(l)), las=1)
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("plot", plot.args)
  if (plot.args$ylab=="") {
    ylab <- if (norm) expression("||"*hat(beta)*"||") else expression(hat(beta))
    mtext(ylab, 2, 3.5, las=1, adj=0)
  }
  abline(h=0,lwd=0.5,col="gray")

  cols <- hcl(h=seq(15,375,len=max(4,n.g+1)),l=60,c=150,alpha=alpha)
  cols <- if (n.g==2) cols[c(1,3)] else cols[1:n.g]
  line.args <- list(col=cols, lwd=1+2*exp(-p/20), lty=1, pch="")
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- l
  line.args$y <- t(Y)
  line.args$col <- line.args$col[g]
  line.args$lty <- rep(line.args$lty, length.out=max(g))
  line.args$lty <- line.args$lty[g]
  do.call("matlines",line.args)
  
  if(!missing(legend.loc)) {
    legend.args <- list(col=cols, lwd=line.args$lwd, lty=line.args$lty, legend=names(x$group.multiplier))
    if (length(new.args)) {
      new.legend.args <- new.args[names(new.args) %in% names(formals(legend))]
      legend.args[names(new.legend.args)] <- new.legend.args
    }
    legend.args$x <- legend.loc
    do.call("legend",legend.args)
  }
}
